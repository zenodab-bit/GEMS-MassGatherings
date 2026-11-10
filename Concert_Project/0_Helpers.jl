# ===========================================================================
# 0_Helpers.jl
#
# Shared foundational code, included first. Sections in this file:
#   1. Numeric distribution helpers
#   2. Age-group helpers
#   3. Event configuration model (Section, Category, EventConfig, Event)
#   4. Config validation (validate_config)
#   5. Result types & aggregation (per-event compartment counts)
#   6. Output formatting & plotting helpers
# ===========================================================================


## === 1. Numeric distribution helpers ===

# (r, p) params for NegativeBinomial. Requires variance > mean.
# use sample_n_contacts, not this directly, for inputs that might violate that.
function negbin_params(mean, std)
    variance = std^2
    p = mean / variance
    r = mean^2 / (variance - mean)
    return r, p
end

# Draws a contact count from (mean, std); falls back to Poisson(mean) when
# variance <= mean, since NegativeBinomial isn't valid there.
function sample_n_contacts(rng, mean, std)
    mean <= 0 && return 0
    variance = std^2
    if variance <= mean
        return rand(rng, Poisson(mean))
    end
    r, p = negbin_params(mean, std)
    return rand(rng, NegativeBinomial(r, p))
end

# (α, β) params for Gamma. Valid for any mean, std > 0, but Gamma is
# unbounded above — clamp to [0,1] if the result is used as a probability
# (see prepare_population in 3_Population.jl).
function gamma_params(mean, std)
    variance = std^2
    β = variance / mean
    α = mean / β
    return α, β
end

# (α, β) params for Beta, which always stays within [0, 1]. Needs
# 0 < mean < 1 and std^2 < mean * (1 - mean); errors otherwise.
function beta_params(mean, std)
    0 < mean < 1 || error("beta_params: mean ($mean) must be between 0 and 1.")
    variance = std^2
    variance < mean * (1 - mean) ||
        error("beta_params: std ($std) too large for mean $mean; a Beta needs std^2 < mean * (1 - mean) = $(mean * (1 - mean)).")
    k = mean * (1 - mean) / variance - 1
    return mean * k, (1 - mean) * k
end


## === 2. Age-group helpers ===

# Classifies a raw age into a group index (1..length(age_boundaries)+1) via
# binary search. age_boundaries must be sorted ascending.
age_group_idx(age, age_boundaries) = searchsortedfirst(age_boundaries, age)

# Formats a group index as a label, e.g. index 1 of [45,65] → "<=45".
function age_group_label_from_idx(idx::Int, age_boundaries)
    idx > length(age_boundaries) && return ">$(age_boundaries[end])"
    return "<=$(age_boundaries[idx])"
end

# Classify and format a raw age in one call.
age_group_label(age, age_boundaries) = age_group_label_from_idx(age_group_idx(age, age_boundaries), age_boundaries)


## === 3. Event configuration model ===

# One sub-group within a draw of a Category (e.g. "seated" vs "standing").
# Sections of the same draw share one date but have independent attendees.
@with_kw struct Section
    id::Int
    name::String = ""
    n_range::Tuple{Int,Int}      # event size drawn uniformly from this range
    mean_contacts::Float64       # mean within-event contacts per attendee
    std_contacts::Float64 = 0.0  # needs std^2 > mean_contacts for NegBin; otherwise falls back to Poisson
end

# A repeatable kind of gathering (e.g. "festival"). Recurs n_draws times,
# each on its own date, each with an independently-sampled attendee list
# per section.
@with_kw struct Category
    id::Int                                     # unique across all categories
    name::String = ""                           # should be unique across categories (used in filenames)
    date_range::Tuple{Int,Int}                  # (first_day, last_day) draws may land on
    sections::Vector{Section}
    n_draws::Int                                # also the max times one person can attend this category
    min_age::Int = 0                            # eligibility filter
    max_age::Int = 999
    age_weights::Vector{Float64} = Float64[]    # overrides EventConfig.age_dist if non-empty
    sex_weights::Vector{Float64} = Float64[]    # overrides EventConfig.sex_dist if non-empty
    core::Float64 = 0.0                         # fraction (0-1) of the smallest draw (all sections together) that attends every draw
    loyalty::Float64 = 0.0                      # guaranteed, demographically-weighted fraction of non-core spots for past non-core attendees; pool never expires, no extra weight for repeat attendance
    cross_section_mean_contacts::Float64 = 0.0  # contacts with attendees of other sections of the same draw; 0 keeps sections isolated
    cross_section_std_contacts::Float64 = 0.0
    # Where core and loyal members go when there are several sections (no
    # effect with a single section):
    #   :event   - placed in a random section at every draw, in proportion to section size
    #   :section - core keep a fixed home section; loyal attendees return to their last section
    section_mode::Symbol = :event
end

# The whole simulation's event setup.
@with_kw struct EventConfig
    categories::Vector{Category}            # ids must be unique
    age_boundaries::Vector{Int}              # sorted ascending, e.g. [45, 65]
    age_dist::Vector{Float64}                # population-wide target age distribution; should sum to 1.0
    sex_dist::Vector{Vector{Float64}}        # population-wide target sex distribution per age group
end

# One concrete gathering: a Section of a draw of a Category, on one date.
# Produced by sample_events() (3_Population.jl) — not user-constructed.
@with_kw struct Event
    id::String                # unique, "$(category_id)_$(draw)_$(section_id)"
    name::String               # human-readable, NOT guaranteed unique
    category_id::Int
    draw_id::Int                # which recurrence (1..n_draws, in date order), not a global counter
    section_id::Int
    date::Int
    n::Int                      # target attendee count
    mean_contacts::Float64
    std_contacts::Float64
    cross_section_mean_contacts::Float64 = 0.0
    cross_section_std_contacts::Float64 = 0.0
end


## === 4. Config validation ===

# Checks an EventConfig for internal consistency — call once, right after
# construction (see end of 1_UserConfig.jl).
# Verifies age_dist/sex_dist have the right length for age_boundaries, age_dist
# sums to ~1.0, each category's age_weights/sex_weights (if given) have
# the right length, min_age <= max_age, date_range and each section's
# n_range aren't inverted, category.core is within [0, 1], n_draws isn't
# negative (0 only warns: it's a valid way to switch a category off), and
# category ids, and section ids within each category, are unique.
function validate_config(event_config::EventConfig)
    n_groups = length(event_config.age_boundaries) + 1

    # duplicate ids would silently merge events (event ids are built from them)
    category_ids = [c.id for c in event_config.categories]
    allunique(category_ids) ||
        error("Category ids must be unique, got $(category_ids).")

    length(event_config.age_dist) == n_groups ||
        error("age_dist has $(length(event_config.age_dist)) entries, expected $n_groups (length(age_boundaries)+1).")

    isapprox(sum(event_config.age_dist), 1.0; atol=1e-6) ||
        @warn "age_dist sums to $(sum(event_config.age_dist)), not 1.0 — event demographics will be silently skewed."

    length(event_config.sex_dist) == n_groups ||
        error("sex_dist has $(length(event_config.sex_dist)) entries, expected $n_groups (one per age group).")

    for (i, sd) in enumerate(event_config.sex_dist)
        length(sd) == 2 || error("sex_dist[$i] has $(length(sd)) entries, expected 2.")
    end

    for category in event_config.categories
        category.n_draws >= 0 ||
            error("Category $(category.name): n_draws ($(category.n_draws)) must not be negative.")
        category.n_draws == 0 &&
            @warn "Category $(category.name): n_draws is 0, so it creates no events."

        section_ids = [s.id for s in category.sections]
        allunique(section_ids) ||
            error("Category $(category.name): section ids must be unique, got $(section_ids).")

        if !isempty(category.age_weights)
            length(category.age_weights) == n_groups ||
                error("Category $(category.name): age_weights has $(length(category.age_weights)) entries, expected $n_groups.")
        end
        if !isempty(category.sex_weights)
            length(category.sex_weights) == 2 ||
                error("Category $(category.name): sex_weights has $(length(category.sex_weights)) entries, expected 2.")
        end
        category.min_age <= category.max_age ||
            error("Category $(category.name): min_age ($(category.min_age)) exceeds max_age ($(category.max_age)).")
        category.date_range[1] <= category.date_range[2] ||
            error("Category $(category.name): date_range $(category.date_range) is invalid (start > end).")
        for section in category.sections
            section.n_range[1] <= section.n_range[2] ||
                error("Category $(category.name), section $(section.id): n_range $(section.n_range) is invalid (start > end).")
        end

        # core > 1 silently oversizes the core group rather than erroring elsewhere
        0.0 <= category.core <= 1.0 ||
            error("Category $(category.name): core ($(category.core)) must be between 0 and 1.")

        category.section_mode in (:event, :section) ||
            error("Category $(category.name): section_mode must be :event or :section, got :$(category.section_mode).")
    end

    println("Config validated OK.")
end


## === 5. Result types & aggregation ===

# Concrete types (not Vector{Any}/Dict{Symbol,NamedTuple}) for per-event
# compartment counts and their aggregated stats — for type stability.
const EventCounts = @NamedTuple{susceptible::Int, infectious::Int, exposed::Int, recovered::Int, absent::Int, same_day_other::Int, infected_at_event::Int}
const EventStat = @NamedTuple{mean::Float64, std::Float64, min::Float64, p25::Float64, median::Float64, p75::Float64, max::Float64}

# Whether someone is dead, severely ill or in hospital on day t, from their
# infection-log dates (exposure = the tick they were infected). Uses the
# same rules as GEMS's is_dead / is_severe / is_hospitalized, which
# absent_from_events (4_Contacts.jl) uses during the simulation: such people
# don't attend events that day.
function absent_on(t, exposure, recovery, death, sev_onset, sev_offset, hosp_admission, hosp_discharge)
    0 <= death <= t && return true
    infected = exposure >= 0 && exposure <= t < max(recovery, death)
    infected || return false
    return (0 <= sev_onset <= t < sev_offset) || (0 <= hosp_admission <= t < hosp_discharge)
end

absent_on(t, row) = absent_on(t, row.tick, row.recovery, row.death, row.severeness_onset,
    row.severeness_offset, row.hospital_admission, row.hospital_discharge)

# Classifies every assigned attendee of every event into one compartment as
# of that event's date, by scanning the infection log once per event.
# Absent (dead, severe or hospitalised that day) comes first; the other
# compartments only count attendees who were actually there.
function analyze_event_population(inf_log::DataFrame, people::DataFrame, events::Vector{Event}, attendees::Dict{String,Vector{Int}})
    results = Dict{String,EventCounts}()

    for event in events
        # get individual ids for this event from people DataFrame
        section_ids = Set{Int32}(people.id[idx] for idx in attendees[event.id])

        # categorize people from infection log
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
                    push!(event_infected_ids, row.id_b)
                end
            end
        end

        # compute counts
        not_susceptible_ids = union(infected_before_ids, same_day_other_ids)
        susceptible = length(setdiff(section_ids, not_susceptible_ids))
        infectious = length(intersect(currently_infectious_ids, section_ids))
        exposed = length(intersect(exposed_before_ids, section_ids))
        recovered = length(intersect(recovered_ids, section_ids))
        absent = length(intersect(absent_ids, section_ids))
        same_day_other = length(intersect(same_day_other_ids, section_ids))
        infected_at_event = length(intersect(event_infected_ids, section_ids))

        event_data::EventCounts = (
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


# Aggregates per-run EventCounts into mean/std/min/median/max/quartiles per
# event per field, across replicates. std is 0.0 with a single replicate.
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


# Per-run attack rate among susceptible attendees of one event, in %:
# infected_at_event / susceptible. Runs with no susceptible attendees are
# skipped (the rate is undefined there), so the result can have fewer
# entries than event_results, or none.
function susceptible_attack_rates(event_results::Vector{Dict{String,EventCounts}}, event_id::String)
    rates = Float64[]
    for r in event_results
        c = r[event_id]
        c.susceptible > 0 && push!(rates, 100 * c.infected_at_event / c.susceptible)
    end
    return rates
end


# Dispersion k of onward infections (offspring) for one run: lower k means
# more superspreading. Counts how many people each infected person went on
# to infect (0 for those who infected nobody) and fits k by method of
# moments, k = mean^2 / (variance - mean); Inf means no overdispersion.
#
# People infected too close to the end are left out, since the simulation
# ends before they finish infecting others. The cutoff is measured from this
# run's own infection log, so it follows whatever the TOML's pathogen
# settings are: for every infected person, the days from infection until
# they can no longer infect anyone, i.e. the earliest of recovery, death or
# hospital admission (hospitalised people can't infect; severe people still
# infect their household until recovery). The cutoff is the `coverage`
# quantile (default 95%) of those durations, rounded up to whole days.
#
# The count pools the whole epidemic, so it also reflects susceptibles
# running out over time.
function offspring_dispersion(inf_log::DataFrame, last_tick::Integer; coverage::Real=0.95)
    offspring = Dict{Int32,Int}()
    for a in inf_log.id_a
        a >= 0 || continue   # seed infections have no infecter (negative id)
        offspring[a] = get(offspring, a, 0) + 1
    end

    # days each person stays able to infect (-1 in the log = never happens)
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


## === 6. Output formatting & plotting helpers ===

# Colours for every plot: the Okabe-Ito palette, chosen to stay
# distinguishable for the common forms of colour blindness. Change a colour
# here to change it everywhere. 6_Analysis.jl makes PLOT_PALETTE the default
# palette, so plots that don't set colours themselves (e.g. GEMS's epidemic
# overview) use it in this order. Black is left out of the palette, since it
# hides the median lines in boxplots and clashes with the event-day lines.
const PLOT_COLORS = (
    blue       = colorant"#0072B2",
    orange     = colorant"#E69F00",
    green      = colorant"#009E73",
    vermillion = colorant"#D55E00",
    purple     = colorant"#CC79A7",
    skyblue    = colorant"#56B4E9",
    yellow     = colorant"#F0E442",   # low contrast on white: last in the palette
    black      = colorant"#000000",
)
const PLOT_PALETTE = [PLOT_COLORS.blue, PLOT_COLORS.orange, PLOT_COLORS.green, PLOT_COLORS.vermillion,
                      PLOT_COLORS.purple, PLOT_COLORS.skyblue, PLOT_COLORS.yellow]

# Saves what's needed to rerun a run into <run_folder>/run_info:
#   - copies of `config_files` (the files the run was configured with)
#   - run_info.toml: the `settings` passed in (seeds, replicates, population
#     files, ...), Julia version and thread count, and the code version
#     (git commit and branch, and whether there were uncommitted changes)
#   - uncommitted_changes.diff: those changes, if any (new files that were
#     never added to git are not included)
# The population and settings files are recorded by path, not copied.
# See README, "Rerunning a run".
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
        "environment" => Dict("julia_version" => string(VERSION), "threads" => Threads.nthreads()),
        "code" => Dict(
            "git_commit" => git_out(`git -C $repo rev-parse HEAD`),
            "git_branch" => git_out(`git -C $repo rev-parse --abbrev-ref HEAD`),
            "uncommitted_changes" => has_changes,
        ),
    )
    open(io -> TOML.print(io, info), joinpath(dir, "run_info.toml"), "w")
    println("Saved: run_info (config copies, seeds, code version)")
end

# Saves the per-event results before they're summarised, so they can be
# reused without rerunning:
#   - events.csv: one row per event (id, category, draw, section, date, size,
#     contact rates)
#   - event_results.csv: one row per replicate and event, with every
#     compartment count (susceptible ... infected_at_event)
#   - attendance.jld2: who was assigned to each event ("attendees": event id
#     => person ids) and each category's core group ("core_groups":
#     category id => person ids). Person ids are the population's `id`.
function save_raw_results(run_folder, events, event_results, attendees, core_groups, people)
    events_df = DataFrame([NamedTuple{fieldnames(Event)}(getfield.(Ref(e), fieldnames(Event))) for e in events])
    CSV.write(joinpath(run_folder, "events.csv"), events_df)

    results_df = DataFrame([merge((replicate=r, event_id=e.id), res[e.id])
                            for (r, res) in enumerate(event_results) for e in events])
    CSV.write(joinpath(run_folder, "event_results.csv"), results_df)

    JLD2.save(joinpath(run_folder, "attendance.jld2"),
        "attendees", Dict(id => people.id[idxs] for (id, idxs) in attendees),
        "core_groups", Dict(cid => people.id[idxs] for (cid, idxs) in core_groups))
    println("Saved: events.csv, event_results.csv, attendance.jld2")
end

# Rounds a number for consistent display in printed tables.
function fmt(x; digits=2)
    round(x, digits=digits)
end

# Column-header row for the per-event metrics table (multievent_analysis.txt).
function print_table_header(io)
    println(io, rpad("Metric", 20), " | ",
        rpad("Mean", 9), " | ", rpad("Std", 9), " | ",
        rpad("Min", 9), " | ", rpad("P25", 9), " | ",
        rpad("Median", 9), " | ", rpad("P75", 9), " | ", "Max")
    println(io, "-"^100)
end

# One row of the per-event metrics table.
function print_metric_row(io, label, m)
    println(io, rpad(label, 20), " | ",
        rpad(fmt(m.mean), 9), " | ", rpad(fmt(m.std), 9), " | ",
        rpad(fmt(m.min), 9), " | ", rpad(fmt(m.p25), 9), " | ",
        rpad(fmt(m.median), 9), " | ", rpad(fmt(m.p75), 9), " | ",
        fmt(m.max))
end

# One dashed vline per unique event date, labeled with every event.id
# sharing that date (avoids duplicate overlapping lines/legend entries).
function add_event_vlines!(p, events)
    seen_dates = Set{Int}()
    for event in events
        event.date in seen_dates && continue
        push!(seen_dates, event.date)

        same_date_ids = [e.id for e in events if e.date == event.date]
        label = length(same_date_ids) > 1 ?
                "Events $(join(same_date_ids, ", ")) (day $(event.date))" :
                "Event $(event.id) (day $(event.date))"

        vline!(p, [event.date],
            linestyle=:dash,
            color=:black,
            label=label,
            linewidth=1.5)
    end
end


println("End Helpers")