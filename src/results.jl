# ===========================================================================
# results.jl
#
# Per-event compartment counts and their summary over replicates, attack
# rates, offspring dispersion, and saving run_info/ and the raw results.
# ===========================================================================

## === Result types and aggregation ===

# Concrete types, for type stability. assigned: the people actually
# assigned in that replicate (below n only if the event couldn't be
# filled); the compartments susceptible to same_day_other add up to it.
const EventCounts = @NamedTuple{assigned::Int, susceptible::Int, infectious::Int, exposed::Int, recovered::Int, absent::Int, same_day_other::Int, infected_at_event::Int}
const EventStat = @NamedTuple{mean::Float64, std::Float64, min::Float64, p25::Float64, median::Float64, p75::Float64, max::Float64}

# Whether someone is dead, severely ill or in hospital on day t, from
# their infection-log dates: the same rules as GEMS's is_dead / is_severe /
# is_hospitalized, which absent_from_events (contacts.jl) uses.
function absent_on(t, exposure, recovery, death, sev_onset, sev_offset, hosp_admission, hosp_discharge)
    0 <= death <= t && return true
    infected = exposure >= 0 && exposure <= t < max(recovery, death)
    infected || return false
    return (0 <= sev_onset <= t < sev_offset) || (0 <= hosp_admission <= t < hosp_discharge)
end

absent_on(t, row) = absent_on(t, row.tick, row.recovery, row.death, row.severeness_onset,
    row.severeness_offset, row.hospital_admission, row.hospital_discharge)

# Classifies every attendee as of the event day, from the infection log.
# `passes`: (person id, day) => stage for infections in later passes. At a
# stage-s venue, earlier stages' infections count as same_day_other.
function analyze_event_population(inf_log::DataFrame, people::DataFrame, events::Vector{Event},
                                  attendees::Dict{String,Vector{Int}};
                                  passes::AbstractDict = Dict{Tuple{Int32,Int16},Int}())
    results = Dict{String,EventCounts}()

    for event in events
        this_pass = max(event.stage_id, 1)
        section_ids = Set{Int32}(people.id[idx] for idx in attendees[event.id])

        infected_before_ids = Set{Int32}()
        currently_infectious_ids = Set{Int32}()
        exposed_before_ids = Set{Int32}()
        recovered_ids = Set{Int32}()
        absent_ids = Set{Int32}()
        same_day_other_ids = Set{Int32}()
        event_infected_ids = Set{Int32}()

        for row in eachrow(inf_log)
            if row.tick < event.date
                push!(infected_before_ids, row.id_b)
                if absent_on(event.date, row)   # includes everyone dead by then
                    push!(absent_ids, row.id_b)
                elseif row.infectiousness_onset <= event.date &&
                   (row.recovery > event.date || row.recovery == -1)
                    push!(currently_infectious_ids, row.id_b)
                elseif row.recovery != -1 && row.recovery <= event.date
                    push!(recovered_ids, row.id_b)
                else
                    push!(exposed_before_ids, row.id_b)
                end
            elseif row.tick == event.date
                if row.setting_type != 'g'
                    push!(same_day_other_ids, row.id_b)
                else
                    row_pass = get(passes, (Int32(row.id_b), Int16(row.tick)), 1)
                    if row_pass == this_pass
                        push!(event_infected_ids, row.id_b)
                    elseif row_pass < this_pass
                        push!(same_day_other_ids, row.id_b)   # at an earlier stage that day
                    end                                       # (a later stage: still susceptible here)
                end
            end
        end

        not_susceptible_ids = union(infected_before_ids, same_day_other_ids)
        susceptible = length(setdiff(section_ids, not_susceptible_ids))
        infectious = length(intersect(currently_infectious_ids, section_ids))
        exposed = length(intersect(exposed_before_ids, section_ids))
        recovered = length(intersect(recovered_ids, section_ids))
        absent = length(intersect(absent_ids, section_ids))
        same_day_other = length(intersect(same_day_other_ids, section_ids))
        infected_at_event = length(intersect(event_infected_ids, section_ids))

        event_data::EventCounts = (
            assigned=length(attendees[event.id]),
            susceptible=susceptible,
            infectious=infectious,
            exposed=exposed,
            recovered=recovered,
            absent=absent,
            same_day_other=same_day_other,
            infected_at_event=infected_at_event
        )

        results[event.id] = event_data
    end

    return results
end


# Mean, std (0.0 for one replicate), min, quartiles and max of every field,
# per event, across replicates.
function aggregate_event_results(results_vector::Vector{Dict{String,EventCounts}})
    aggregated = Dict{String,Dict{Symbol,EventStat}}()
    first_result = results_vector[1]

    for (event_id, _) in first_result
        fields = keys(results_vector[1][event_id])
        aggregated[event_id] = Dict{Symbol,EventStat}()
        for field in fields
            values = Float64[getfield(results_vector[r][event_id], field)
                             for r in 1:length(results_vector)]
            stat::EventStat = (
                mean=mean(values),
                std=length(values) > 1 ? std(values) : 0.0,
                min=minimum(values),
                p25=quantile(values, 0.25),
                median=median(values),
                p75=quantile(values, 0.75),
                max=maximum(values)
            )
            aggregated[event_id][field] = stat
        end
    end

    return aggregated
end


# Per run, infected_at_event / susceptible in %, for one event. Runs with no
# susceptible attendees are skipped, so there may be fewer rates than runs.
function susceptible_attack_rates(event_results::Vector{Dict{String,EventCounts}}, event_id::String)
    rates = Float64[]
    for r in event_results
        c = r[event_id]
        c.susceptible > 0 && push!(rates, 100 * c.infected_at_event / c.susceptible)
    end
    return rates
end


# Dispersion k of onward infections (method of moments; lower k, more
# superspreading). Leaves out people infected in the last `cutoff` days:
# the `coverage` quantile of how long people can infect, from this run.
function offspring_dispersion(inf_log::DataFrame, last_tick::Integer; coverage::Real=0.95)
    offspring = Dict{Int32,Int}()
    for a in inf_log.id_a
        a >= 0 || continue   # seed infections have no infecter (negative id)
        offspring[a] = get(offspring, a, 0) + 1
    end

    # days each person can infect (-1 in the log = never happens)
    windows = Int[]
    tick, rec, dth, hosp = inf_log.tick, inf_log.recovery, inf_log.death, inf_log.hospital_admission
    for i in eachindex(tick)
        stop = typemax(Int)
        rec[i] >= 0 && (stop = min(stop, rec[i]))
        dth[i] >= 0 && (stop = min(stop, dth[i]))
        hosp[i] >= 0 && (stop = min(stop, hosp[i]))
        stop == typemax(Int) && continue
        push!(windows, stop - tick[i])
    end
    isempty(windows) && return (k=NaN, mean=NaN, n=0, cutoff=0)
    cutoff = ceil(Int, quantile(windows, coverage))

    included = tick .<= last_tick - cutoff
    counts = [get(offspring, b, 0) for b in inf_log.id_b[included]]
    length(counts) < 2 && return (k=NaN, mean=NaN, n=length(counts), cutoff=cutoff)

    m = mean(counts)
    v = var(counts)
    k = v > m ? m^2 / (v - m) : Inf
    return (k=k, mean=m, n=length(counts), cutoff=cutoff)
end


# Saves what's needed to rerun a run into run_info/: config copies,
# run_info.toml and uncommitted changes (see docs/reproducibility.md).
function save_run_info(run_folder, config_files::Vector{String}; settings...)
    dir = joinpath(run_folder, "run_info")
    mkpath(dir)
    for f in config_files
        cp(f, joinpath(dir, basename(f)); force=true)
    end

    repo = normpath(joinpath(@__DIR__, ".."))
    git_out(cmd) = try
        strip(read(pipeline(cmd; stderr=devnull), String))
    catch
        "unknown (git not available)"
    end
    diff = git_out(`git -C $repo diff HEAD`)
    has_changes = !isempty(diff) && !startswith(diff, "unknown")
    has_changes && write(joinpath(dir, "uncommitted_changes.diff"), diff * "\n")

    info = Dict(
        "saved_at" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "run" => Dict(string(k) => (isnothing(v) ? "nothing" : v) for (k, v) in settings),
        # interactive threads matter too: plain `julia` (1.12+) starts one,
        # `julia -t N` none, and GEMS's random numbers differ between them
        "environment" => Dict("julia_version" => string(VERSION), "threads" => Threads.nthreads(),
                              "interactive_threads" => Threads.nthreads(:interactive)),
        "code" => Dict(
            "git_commit" => git_out(`git -C $repo rev-parse HEAD`),
            "git_branch" => git_out(`git -C $repo rev-parse --abbrev-ref HEAD`),
            "uncommitted_changes" => has_changes,
        ),
    )
    open(io -> TOML.print(io, info), joinpath(dir, "run_info.toml"), "w")
    println("Saved: run_info (config copies, seeds, code version)")
end

# Saves events.csv, event_results.csv (with history and future columns
# from `labels`) and attendance.jld2 (person ids, per replicate if they vary).
function save_raw_results(run_folder, events, event_results, setups, people; labels = nothing)
    events_df = DataFrame([NamedTuple{fieldnames(Event)}(getfield.(Ref(e), fieldnames(Event))) for e in events])
    # the stage column only matters with staged categories
    all(==(0), events_df.stage_id) && select!(events_df, Not(:stage_id))
    CSV.write(joinpath(run_folder, "events.csv"), events_df)

    label(r) = isnothing(labels) ? (;) : labels[r]
    results_df = DataFrame([merge((replicate=r,), label(r), (event_id=e.id,), res[e.id])
                            for (r, res) in enumerate(event_results) for e in events])
    CSV.write(joinpath(run_folder, "event_results.csv"), results_df)

    with_ids(s) = (Dict(id => people.id[idxs] for (id, idxs) in s.attendees),
                   Dict(cid => people.id[idxs] for (cid, idxs) in s.core_groups))
    if length(setups) == 1
        attendees, core_groups = with_ids(setups[1])
    else
        per_replicate = with_ids.(setups)
        attendees, core_groups = first.(per_replicate), last.(per_replicate)
    end
    JLD2.save(joinpath(run_folder, "attendance.jld2"), "attendees", attendees, "core_groups", core_groups)
    println("Saved: events.csv, event_results.csv, attendance.jld2")
end

