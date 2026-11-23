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
#   7. Running a scenario (setup_streams, run_scenario, run_comparison)
#   8. Sweeps (with_date, with_date_range, run_sweep)
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

# Common supertype of the two kinds of category (Category and
# StagedCategory below), so an EventConfig can hold both.
abstract type AbstractCategory end

# A repeatable kind of gathering (e.g. "festival"). Recurs n_draws times,
# each on its own date, each with an independently-sampled attendee list
# per section.
@with_kw struct Category <: AbstractCategory
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

# One place within a stage of a StagedCategory (e.g. one restaurant, or one
# stadium): a single group, where the attendees present meet each other,
# with no contacts to other venues.
@with_kw struct Venue
    name::String = ""
    n_range::Tuple{Int,Int}      # size drawn uniformly from this range, anew at every draw
    mean_contacts::Float64       # mean contacts per attendee within the venue
    std_contacts::Float64 = 0.0  # needs std^2 > mean_contacts for NegBin; otherwise falls back to Poisson
end

# One stage of a StagedCategory (e.g. "restaurants" or "match"): its venues
# are on the same day as the other stages of the draw. The stage's
# attendees are drawn together, then spread over its venues in proportion
# to their sizes, anew at every draw.
@with_kw struct Stage
    name::String = ""
    venues::Vector{Venue}
    overlap::Float64 = 0.0                      # share (0-1) of this stage's attendees taken from the previous stage's attendees of the same draw; must be 0 for the first stage
    min_age::Int = 0                            # eligibility filter
    max_age::Int = 999
    age_weights::Vector{Float64} = Float64[]    # overrides EventConfig.age_dist if non-empty
    sex_weights::Vector{Float64} = Float64[]    # overrides EventConfig.sex_dist if non-empty
end

# A gathering made of stages on the same day (e.g. restaurants, then a
# football match), where someone can attend several stages of a draw:
# the overlap of a stage says how many come from the previous one; the rest
# are people at no earlier stage of the draw. Recurs n_draws times, each on
# its own date. No core or loyalty: each draw draws its attendees afresh.
@with_kw struct StagedCategory <: AbstractCategory
    id::Int                                     # unique across all categories
    name::String = ""                           # should be unique across categories (used in filenames)
    date_range::Tuple{Int,Int}                  # (first_day, last_day) draws may land on
    stages::Vector{Stage}                       # in order: overlap is always with the previous one
    n_draws::Int
end

# The whole simulation's event setup.
@with_kw struct EventConfig
    categories::Vector{AbstractCategory}    # Category and/or StagedCategory; ids must be unique
    age_boundaries::Vector{Int}              # sorted ascending, e.g. [45, 65]
    age_dist::Vector{Float64}                # population-wide target age distribution; should sum to 1.0
    sex_dist::Vector{Vector{Float64}}        # population-wide target sex distribution per age group
end

# One concrete gathering, on one date: a Section of a draw of a Category,
# or a Venue of a stage of a draw of a StagedCategory.
# Produced by sample_events() (3_Population.jl) — not user-constructed.
@with_kw struct Event
    id::String                # unique: "$(category_id)_$(draw)_$(section_id)", or for a staged category "$(category_id)_$(draw)_$(stage_id)_$(venue)"
    name::String               # human-readable, NOT guaranteed unique
    category_id::Int
    draw_id::Int                # which recurrence (1..n_draws, in date order), not a global counter
    section_id::Int             # the section; for a staged category, the venue's number within its stage
    date::Int
    n::Int                      # target attendee count
    mean_contacts::Float64
    std_contacts::Float64
    cross_section_mean_contacts::Float64 = 0.0
    cross_section_std_contacts::Float64 = 0.0
    stage_id::Int = 0           # the stage (1, 2, ...) of a staged category; 0 for a normal Category
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
        category.date_range[1] <= category.date_range[2] ||
            error("Category $(category.name): date_range $(category.date_range) is invalid (start > end).")

        if category isa StagedCategory
            validate_stages(category, n_groups)
            continue
        end

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

# The checks for a StagedCategory's stages and venues: at least one stage,
# each with at least one venue; the first stage's overlap is 0 (there's no
# stage before it) and every overlap is within [0, 1]; weights of the right
# length, min_age <= max_age, and venue n_ranges not inverted.
function validate_stages(category::StagedCategory, n_groups::Integer)
    label = "Staged category $(category.name)"
    isempty(category.stages) && error("$label: needs at least one stage.")
    for (s, stage) in enumerate(category.stages)
        where = "$label, stage $s ($(stage.name))"
        isempty(stage.venues) && error("$where: needs at least one venue.")
        0.0 <= stage.overlap <= 1.0 || error("$where: overlap ($(stage.overlap)) must be between 0 and 1.")
        s == 1 && stage.overlap > 0 &&
            error("$where: the first stage has no stage before it, so its overlap must be 0, got $(stage.overlap).")
        isempty(stage.age_weights) || length(stage.age_weights) == n_groups ||
            error("$where: age_weights has $(length(stage.age_weights)) entries, expected $n_groups.")
        isempty(stage.sex_weights) || length(stage.sex_weights) == 2 ||
            error("$where: sex_weights has $(length(stage.sex_weights)) entries, expected 2.")
        stage.min_age <= stage.max_age ||
            error("$where: min_age ($(stage.min_age)) exceeds max_age ($(stage.max_age)).")
        for (v, venue) in enumerate(stage.venues)
            venue.n_range[1] <= venue.n_range[2] ||
                error("$where, venue $v ($(venue.name)): n_range $(venue.n_range) is invalid (start > end).")
        end
    end
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
#
# Staged categories: `passes` is the run's record of infections in later
# GlobalSetting passes ((person id, day) => stage, from run_scenario's
# stage_passes); event infections not in it happened in pass 1. For a venue
# of stage s, only infections in stage s count as infected at the event;
# people infected at an earlier stage that day count with same_day_other
# (they're no longer susceptible when stage s comes), and people infected at
# a later stage were still susceptible here. Normal events are pass 1.
function analyze_event_population(inf_log::DataFrame, people::DataFrame, events::Vector{Event},
                                  attendees::Dict{String,Vector{Int}};
                                  passes::AbstractDict = Dict{Tuple{Int32,Int16},Int}())
    results = Dict{String,EventCounts}()

    for event in events
        this_pass = max(event.stage_id, 1)
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
                    row_pass = get(passes, (Int32(row.id_b), Int16(row.tick)), 1)
                    if row_pass == this_pass
                        push!(event_infected_ids, row.id_b)
                    elseif row_pass < this_pass
                        push!(same_day_other_ids, row.id_b)   # infected at an earlier stage that day
                    end                                       # (a later stage: still susceptible here)
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
        # interactive threads too: plain `julia` (1.12+) starts one, `julia -t N` none,
        # and GEMS's random numbers differ between the two even with the same `threads`
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

# Saves the per-event results before they're summarised, so they can be
# reused without rerunning:
#   - events.csv: one row per event (id, category, draw, section, date, size,
#     contact rates)
#   - event_results.csv: one row per replicate and event, with every
#     compartment count (susceptible ... infected_at_event)
#   - attendance.jld2: who was assigned to each event ("attendees": event id
#     => person ids) and each category's core group ("core_groups":
#     category id => person ids). Person ids are the population's `id`.
#     When attendance varies per replicate (several setups), both hold one
#     such Dict per replicate instead, in replicate order.
#   With picked histories or branching, `labels` gives each run's history
#   (replicate number) and future, added as columns to event_results.csv.
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


## === 7. Running a scenario ===

# The setup's random streams, all derived from the setup seed (the TOML's
# Simulation seed). Each setup step gets its own stream, so changing one
# step doesn't reshuffle the others: e.g. changing an event's size or
# n_draws leaves everyone's transmission_prob unchanged.
#   events     - event dates and sizes
#   traits     - per-person traits (transmission_prob)
#   assign     - who attends what
#   population - a synthetic population (only test.jl uses it)
# Created fresh on every call, so the same seed always gives the same streams.
function setup_streams(seed::Integer)
    setup_rng = Xoshiro(seed)
    return (events     = Xoshiro(rand(setup_rng, UInt64)),
            traits     = Xoshiro(rand(setup_rng, UInt64)),
            assign     = Xoshiro(rand(setup_rng, UInt64)),
            population = Xoshiro(rand(setup_rng, UInt64)))
end

# The stream replicate r draws its own attendance (kind = :attendance) or
# transmission_prob (kind = :transmission) from, when they vary per
# replicate (vary_per_replicate in run_scenario). Derived from the setup
# seed with draws after setup_streams' four, so those stay unchanged.
function replicate_stream(setup_seed::Integer, kind::Symbol, r::Integer)
    setup_rng = Xoshiro(setup_seed)
    draws = [rand(setup_rng, UInt64) for _ in 1:6]   # 1-4: setup_streams
    base_seed = kind == :transmission ? draws[5] : kind == :attendance ? draws[6] :
                error("Unknown kind $kind: use :attendance or :transmission.")
    return Xoshiro(base_seed + UInt64(r))
end

# Calls f with the prepared population in the form GEMS needs. GEMS only
# loads a settings file when the population is given as a file path, not as
# a table. So with a settings file, the prepared population (with the event
# columns) is written to a temporary file and passed by path, then deleted
# afterwards, even if f fails; reading it back first makes sure no event
# column changed on the way. Without one, the table is passed directly,
# which is faster (no file to write, and none to read per replicate).
# GEMS_COLUMNS (3_Population.jl) are the columns GEMS copies onto its own
# individuals, for 4_Contacts.jl and 5_Transmission.jl to read.
function with_population_for_gems(f, people, settings_file)
    isnothing(settings_file) && return f(people)
    file = tempname() * ".jld2"
    try
        JLD2.save(file, "data", people)
        let saved = JLD2.load(file, "data")
            for col in GEMS_COLUMNS
                isequal(saved[!, col], people[!, col]) ||
                    error("Column $col changed when writing the population to $file.")
            end
        end
        return f(file)
    finally
        rm(file; force=true)
    end
end

# Branching: makes future f of history h continue with its own random
# numbers from day branch_day on. GEMS keeps one random generator per
# thread (sim.rngs); at the end of day branch_day - 1 (GEMS runs the
# stepmod before moving to the next day) they're re-seeded in place, with
# copy!, so anything holding on to one of them sees the new state too.
# Future 1 isn't re-seeded: it's the history's own continuation, the same
# as without branching. The seeds come from the GEMS seed, h, f and the
# thread's slot k, so the same settings always give the same futures
# (Xoshiro hashes its seed, so nearby numbers give unrelated streams).
future_seed(gems_seed, h, f, k) =
    reinterpret(UInt64, Int64(gems_seed)) ⊻ (UInt64(h) << 40) ⊻ (UInt64(f) << 20) ⊻ UInt64(k)

function branch!(sim, branch_day::Integer, gems_seed, h::Integer, f::Integer)
    f == 1 && return sim
    reseed!(s) = for (k, g) in enumerate(s.rngs)
        copy!(g, Xoshiro(future_seed(gems_seed, h, f, k)))
    end
    if branch_day == 0
        reseed!(sim)   # before the first day
    else
        previous = sim.stepmod   # keep the stepmod already set (daily_global_setting!)
        sim.stepmod = s -> (previous(s); tick(s) == branch_day - 1 && reseed!(s); nothing)
    end
    return sim
end

# The base population table from a .jld2 file (the table under "data").
load_population(file::AbstractString) = JLD2.load(file)["data"]

# The run_scenario settings `kwargs` with a population file replaced by its
# table, loaded once. Used by run_comparison and run_sweep, which call
# run_scenario several times: without it, every run would read the file
# again (about 10 s for Saalekreis). Same results either way, since
# run_scenario always works on its own copy of the table.
function population_loaded_once(kwargs)
    settings = values(kwargs)   # a NamedTuple
    haskey(settings, :population) && settings.population isa AbstractString || return settings
    return merge(settings, (; population = load_population(settings.population)))
end

# Runs one scenario from start to finish: samples the events, prepares the
# population and assigns attendees, runs the GEMS batch, then analyses,
# saves and (optionally) validates the results. 2_Interface.jl and test.jl
# both call this, so the steps exist only here.
#
#   event_config      - the event categories (1_UserConfig.jl)
#   transmission_dist - each person's transmission_prob is drawn from it
#   configfile        - the GEMS TOML
#   population        - the base population: a .jld2 file holding a table
#                       under "data", or the table itself. Never changed:
#                       every call works on its own fresh copy.
#   settings_file     - the matching GEMS settings file, or nothing
#   n_simulations     - number of replicates
#   setup_seed        - seed of the setup streams; defaults to the TOML's
#                       Simulation seed
#   gems_seed         - GEMS's seed: nothing for a new random one, or a number
#   run_validation    - whether to run 7_Validation.jl's checks
#   events_on         - false: everything is set up as usual, but no event
#                       takes place (see daily_global_setting!, 4_Contacts.jl)
#   vary_per_replicate - what each replicate draws for itself instead of
#                       sharing: [] (nothing, the default), [:attendance],
#                       [:transmission] (transmission_prob) or both. The
#                       events themselves (dates, sizes) are always shared.
#   histories         - which replicates to run, by number: nothing (all of
#                       1:n_simulations) or e.g. [7], which runs replicate 7
#                       exactly as in a run of 7 or more replicates with the
#                       same seeds and settings
#   branch_day, n_futures - branching: every history runs n_futures times,
#                       identical up to day branch_day - 1 and with its own
#                       random numbers from branch_day on (see branch!).
#                       branch_day = nothing: no branching
#   label             - the GEMS batch label
#   results_name      - folder name inside Results/; by default
#                       run_<n>sims_<date and time>
#   run_info_files    - files copied into run_info/ (see save_run_info)
#   run_info          - extra settings written to run_info.toml
#
# Returns a NamedTuple with everything the run produced: events, people,
# attendees, core_groups (replicate 1's, when they vary), setups (one per
# replicate when something varies, otherwise one shared by all: each with
# attendees, core_groups and transmission_prob), bd (the GEMS BatchData),
# event_results, aggregated, stage_passes (per run: (person id, day) =>
# stage, for infections at later stages of staged categories), run_folder,
# setup_seed, gems_seed_used and last_day.
function run_scenario(event_config::EventConfig;
                      transmission_dist,
                      configfile::AbstractString,
                      population,
                      settings_file = nothing,
                      n_simulations::Integer,
                      setup_seed::Integer = TOML.parsefile(configfile)["Simulation"]["seed"],
                      gems_seed = nothing,
                      run_validation::Bool = true,
                      events_on::Bool = true,
                      vary_per_replicate = Symbol[],
                      branch_day = nothing,
                      n_futures::Integer = 1,
                      histories = nothing,
                      label::AbstractString = "MultiEvent simulation",
                      results_name = nothing,
                      run_info_files::Vector{String} = [configfile],
                      run_info = (;))

    isnothing(gems_seed) || gems_seed isa Integer ||
        error("gems_seed must be nothing or an integer, got $(repr(gems_seed)).")
    streams = setup_streams(setup_seed)

    # --- events ---
    events = sample_events(event_config, streams.events)
    # sorted for readable output; assign_events! orders dates itself
    sort!(events, by = e -> (e.date, e.category_id, e.draw_id, e.section_id))

    println("\n=== Sampled Events ===")
    for e in events
        println("Event $(e.id) ($(e.name)) — date: $(e.date), n: $(e.n), contacts: $(e.mean_contacts)")
    end

    sim_length = TOML.parsefile(configfile)["Simulation"]["StopCriterion"]["parameters"]["limit"]
    last_day = sim_length - 1   # GEMS simulates days 0 to limit - 1

    for e in events
        if e.date > last_day
            @warn "Event $(e.id) on day $(e.date) is after the last simulated day ($last_day; " *
                  "GEMS simulates days 0 to StopCriterion limit - 1). Results will be inaccurate."
        end
    end

    vary = Symbol.(collect(vary_per_replicate))
    all(in((:attendance, :transmission)), vary) && allunique(vary) ||
        error("vary_per_replicate may only contain :attendance and :transmission, got $(repr(vary_per_replicate)).")

    # histories (which replicates) and branching (futures per history)
    isnothing(histories) ||
        (histories isa AbstractVector && !isempty(histories) && all(h -> h isa Integer && h >= 1, histories) &&
         allunique(histories)) ||
        error("histories must be nothing or a list of different replicate numbers (1, 2, ...), got $(repr(histories)).")
    n_futures >= 1 || error("n_futures must be at least 1, got $n_futures.")
    n_futures == 1 || !isnothing(branch_day) ||
        error("n_futures = $n_futures needs a branch_day: the day from which the futures differ.")
    isnothing(branch_day) || (branch_day isa Integer && 0 <= branch_day <= last_day) ||
        error("branch_day must be nothing or a day from 0 to $last_day (the last simulated day), got $(repr(branch_day)).")
    history_list = isnothing(histories) ? collect(1:n_simulations) : collect(Int, histories)
    # one entry per GEMS run, in run order
    run_labels = [(history = h, future = f) for h in history_list for f in 1:n_futures]
    labelled = !isnothing(branch_day) || !isnothing(histories)
    # staged categories need each run's record of which stage infections
    # happened in, so their runs go one by one too
    staged = any(c -> c isa StagedCategory, event_config.categories)

    transmission_func = SettingRate()

    # A GEMS batch of n_runs replicates on the prepared population `pop`
    # (a table, or the file it was written to), whose event columns are in
    # the table `people`.
    make_batch(pop, n_runs, people) = Batch(
        n_runs = n_runs,
        configfile = configfile,
        population = pop,
        settingsfile = settings_file,
        ind_extension = GEMS_COLUMNS,
        transmission_function = transmission_func,
        label = label,
        # GlobalSetting holds only each day's event attendees (see 4_Contacts.jl);
        # total event attendances, checked again inside GEMS by daily_global_setting!
        setup = let expected_attendances = sum(length, people.event_dates)
            sim -> daily_global_setting!(sim, expected_attendances; events_on)
        end
    )

    # Picked here rather than by GEMS so the value is known and can be saved
    # with the results (GEMS doesn't report a seed it chose itself).
    gems_seed_used = isnothing(gems_seed) ? rand(RandomDevice(), 0:typemax(Int64)) : Int64(gems_seed)
    Threads.nthreads() == 1 &&
        println("\nRunning on 1 thread. Starting Julia with more, e.g. `julia -t 8 ...`, makes the " *
                "simulation about twice as fast (see README, \"Speed\").")
    println("\nGEMS seed: $gems_seed_used, $(Threads.nthreads()) thread(s). " *
            "To repeat this run, set gems_seed = $gems_seed_used and use the same number of threads.")

    # --- population, attendance and GEMS ---
    # prepare_population and assign_events! add columns to the table, so each
    # call works on a fresh copy and the caller's table stays as it was
    base = population isa AbstractString ? load_population(population) : copy(population)

    # per run: (person id, day) => stage, for infections in the extra
    # GlobalSetting passes of staged categories (see daily_global_setting!)
    stage_passes = Dict{Tuple{Int32,Int16},Int}[]
    if isempty(vary) && !labelled && !staged
        # one setup, shared by every replicate: one GEMS batch
        people, age_groups, sex_levels = prepare_population(event_config, streams.traits, base, transmission_dist)
        attendees, core_groups = assign_events!(people, events, event_config, streams.assign)
        bd = with_population_for_gems(people, settings_file) do pop
            BatchData(make_batch(pop, n_simulations, people); rd_style = "EssentialResultData", seed = gems_seed_used)
        end
        setups = [(; attendees, core_groups, transmission_prob = people.transmission_prob)]
        append!(stage_passes, [Dict{Tuple{Int32,Int16},Int}() for _ in 1:n_simulations])
    else
        # The runs go one by one, the way GEMS's own batch runs them (the same
        # seed for replicate h as in a normal run), collected into one
        # BatchData. Needed when every replicate draws its own attendance
        # and/or transmission_prob (GEMS then gets each its own population),
        # for branching, and for picking histories. Only the slim part of each
        # setup is kept, to save memory.
        sim_seeds = let rng = Xoshiro(gems_seed_used)
            [GEMS.gems_rand(rng, 0:typemax(Int64)) for _ in 1:maximum(history_list)]
        end
        # the setup of history h: shared by all, or its own (vary_per_replicate)
        function make_setup(h)
            fresh = setup_streams(setup_seed)   # streams from their start, as for replicate-shared setups
            traits_rng = :transmission in vary ? replicate_stream(setup_seed, :transmission, h) : fresh.traits
            assign_rng = :attendance in vary ? replicate_stream(setup_seed, :attendance, h) : fresh.assign
            people_h, _, _ = prepare_population(event_config, traits_rng, copy(base), transmission_dist)
            attendees_h, core_groups_h = assign_events!(people_h, events, event_config, assign_rng)
            return people_h, (; attendees = attendees_h, core_groups = core_groups_h,
                                transmission_prob = people_h.transmission_prob)
        end
        shared = isempty(vary) ? make_setup(1) : nothing
        processor = GEMS.BatchProcessor(; keep_rundata = true, master_seed = gems_seed_used)
        setups = []
        people = nothing
        for h in history_list
            people_h, setup_h = isnothing(shared) ? make_setup(h) : shared
            # the first history's full table is kept for saving and validation (ids are the same in all)
            isnothing(people) && (people = people_h)

            with_population_for_gems(people_h, settings_file) do pop
                for f in 1:n_futures
                    println("\n=== Replicate $h" * (n_futures > 1 ? ", future $f of $n_futures" : "") *
                            (isempty(vary) ? "" : ": its own $(join(string.(vary), " and "))") * " ===")
                    push!(setups, setup_h)
                    b = make_batch(pop, 1, people_h)
                    sim = Simulation(; only(GEMS.simconfigs(b))..., seed = sim_seeds[h])
                    only(GEMS.simsetups(b))(sim)
                    isnothing(branch_day) || branch!(sim, branch_day, gems_seed_used, h, f)
                    run!(sim)
                    push!(stage_passes, copy(settings(sim, GlobalSetting)[1].contact_sampling_method.pass_infections))
                    GEMS.accumulate!(processor, PostProcessor(sim); rd_style = "EssentialResultData")
                end
            end
        end
        bd = BatchData(processor; style = "DefaultBatchData")
        isempty(vary) && (setups = [shared[2]])   # one setup for all runs
        attendees, core_groups = setups[1].attendees, setups[1].core_groups
    end
    # replicate r's setup (the same one for all replicates when nothing varies)
    setup_of(r) = setups[length(setups) == 1 ? 1 : r]

    # --- per-event results ---
    event_results = Vector{Dict{String,EventCounts}}()
    for (r, rd) in enumerate(runs(bd))
        push!(event_results, analyze_event_population(infections(rd), people, events, setup_of(r).attendees;
                                                      passes = stage_passes[r]))
    end
    aggregated = aggregate_event_results(event_results)

    # --- saving ---
    # always the project's Results folder, wherever Julia was started from
    name = something(results_name,
        "run_$(length(run_labels))sims_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))")
    run_folder = normpath(joinpath(@__DIR__, "..", "Results", name))
    mkpath(run_folder)

    # everything needed to rerun this run (see README, "Rerunning a run")
    save_run_info(run_folder, run_info_files;
        n_simulations, run_validation, settings_file, events_on,
        vary_per_replicate = string.(vary), branch_day, n_futures, histories, run_info...,
        setup_seed, gems_seed = gems_seed_used)
    # the per-event results before they're summarised
    save_raw_results(run_folder, events, event_results, setups, people;
                     labels = labelled ? run_labels : nothing)

    println("\nSimulation complete.")
    run_analysis(; aggregated, event_results, events, event_config, bd, run_folder,
                 setup_seed, gems_seed_used, last_day)

    # validation output goes to validation.txt, then is shown on screen too
    if run_validation
        run_validation_checks(; people, events, event_config, attendees,
                              core_groups, aggregated, event_results, bd, run_folder, setups)
    end

    return (; events, people, attendees, core_groups, setups, bd, event_results, aggregated,
              stage_passes, run_folder, setup_seed, gems_seed_used, last_day)
end

# Runs a scenario twice, with and without its events, and compares the two
# (compare_runs, 6_Analysis.jl). Both runs get the same setup seed (same
# events, attendees and transmission_prob) and the same GEMS seed, so
# replicate 1 of one run pairs with replicate 1 of the other, and so on:
# each pair is identical up to the first event day, so what differs after
# it comes from the events.
#
# Takes the same arguments as run_scenario. Results go to
# Results/<results_name>, by default compare_<n>sims_<date and time>, with
# one subfolder per run (with_events/, without_events/) and the comparison
# files next to them. Validation (if switched on) runs on the run with
# events only: its checks assume the events took place.
#
# Returns (; with, without, comparison, folder): the two run_scenario
# results, the per-pair comparison table, and the folder.
function run_comparison(event_config::EventConfig;
                        n_simulations::Integer,
                        gems_seed = nothing,
                        run_validation::Bool = true,
                        results_name = nothing,
                        kwargs...)

    isnothing(gems_seed) || gems_seed isa Integer ||
        error("gems_seed must be nothing or an integer, got $(repr(gems_seed)).")
    # picked once here, so both runs use the same one
    gems_seed_used = isnothing(gems_seed) ? rand(RandomDevice(), 0:typemax(Int64)) : Int64(gems_seed)

    name = something(results_name,
        "compare_$(n_simulations)sims_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))")
    folder = normpath(joinpath(@__DIR__, "..", "Results", name))

    kwargs = population_loaded_once(kwargs)

    println("\n######## Run 1 of 2: with events ########")
    with = run_scenario(event_config; kwargs..., n_simulations, gems_seed = gems_seed_used,
                        run_validation, results_name = joinpath(name, "with_events"))

    println("\n######## Run 2 of 2: without events ########")
    without = run_scenario(event_config; kwargs..., n_simulations, gems_seed = gems_seed_used,
                           run_validation = false, events_on = false,
                           results_name = joinpath(name, "without_events"))

    comparison = compare_runs(with, without, folder)
    return (; with, without, comparison, folder)
end


## === 8. Sweeps ===

# A copy of event_config with category `category_id`'s date_range set to
# `date_range`; everything else unchanged. For sweeps over when events happen.
function with_date_range(event_config::EventConfig, category_id::Integer, date_range::Tuple{Int,Int})
    i = findfirst(c -> c.id == category_id, event_config.categories)
    i === nothing && error("No category with id $category_id in event_config.")
    categories = copy(event_config.categories)
    c = categories[i]
    categories[i] = typeof(c)(c; date_range)   # Category or StagedCategory, the rest unchanged
    return EventConfig(event_config; categories)
end

# The same, with every draw of the category on one day. Only for categories
# with n_draws = 1: draws of one category can't share a day (use
# with_date_range to move the window of a recurring category instead).
function with_date(event_config::EventConfig, category_id::Integer, day::Integer)
    i = findfirst(c -> c.id == category_id, event_config.categories)
    i === nothing && error("No category with id $category_id in event_config.")
    n = event_config.categories[i].n_draws
    n <= 1 || error("Category $category_id has n_draws = $n, and its draws can't all be on day $day. " *
                    "Use with_date_range to move its date window instead.")
    return with_date_range(event_config, category_id, (Int(day), Int(day)))
end

# Runs the comparison with and without events (run_comparison) once per
# value in sweep_values, to see how the events' effect depends on one
# setting. `sweep` says what changes: sweep(value) returns the run_scenario
# settings to use for that value, as a NamedTuple, e.g.
#   sweep(day) = (; event_config = with_date(event_config, 1, day))
#   sweep(m)   = (; transmission_dist = Gamma(gamma_params(m, 0.1)...))
# Everything not returned stays as passed in (kwargs, the same as for
# run_scenario). sweep_label names the setting in the output.
#
# Every value uses the same setup seed and the same GEMS seed, so values
# differ only in the swept setting, not in luck.
#
# Shortcut: if sweep only ever returns event_config, the run without events
# is the same for every value (nobody goes to any event, and the
# population, transmission_prob and seeds don't change), so it runs once
# and every value is compared against it. Otherwise each value gets its own.
#
# Results go to Results/<results_name>, by default
# sweep_<n>sims_<date and time>: one folder per value (value_<v>/, holding
# with_events/, its comparison files and, without the shortcut, its own
# without_events/), the shared without_events/ (with the shortcut), and the
# summary over all values: sweep.txt, sweep.csv, sweep.png.
#
# Returns (; summary, folder): the per-value summary table and the folder.
function run_sweep(event_config::EventConfig;
                   sweep::Function,
                   sweep_values::AbstractVector,
                   sweep_label::AbstractString = "value",
                   n_simulations::Integer,
                   gems_seed = nothing,
                   run_validation::Bool = true,
                   results_name = nothing,
                   kwargs...)

    isempty(sweep_values) && error("sweep_values is empty.")
    allunique(sweep_values) || error("sweep_values has repeated values: $(sweep_values).")
    isnothing(gems_seed) || gems_seed isa Integer ||
        error("gems_seed must be nothing or an integer, got $(repr(gems_seed)).")
    # picked once here, so every run uses the same one
    gems_seed_used = isnothing(gems_seed) ? rand(RandomDevice(), 0:typemax(Int64)) : Int64(gems_seed)

    kwargs = population_loaded_once(kwargs)
    overrides = [sweep(v) for v in sweep_values]
    for (v, o) in zip(sweep_values, overrides)
        o isa NamedTuple ||
            error("sweep($v) returned $(typeof(o)); it must return settings like (; event_config = ...).")
    end
    events_only = all(o -> keys(o) == (:event_config,), overrides)

    name = something(results_name,
        "sweep_$(n_simulations)sims_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))")
    folder = normpath(joinpath(@__DIR__, "..", "Results", name))
    n_runs = length(sweep_values) + (events_only ? 1 : length(sweep_values))
    println("\n=== Sweep: $(length(sweep_values)) values of $sweep_label, $n_runs runs of " *
            "$n_simulations replicate(s) ($(events_only ? "one run without events for all values" :
                                              "one run without events per value")) ===")

    # settings shared by every run; for each value, its own settings win
    shared = (; kwargs..., n_simulations, gems_seed = gems_seed_used)
    settings_for(o) = merge(shared, Base.structdiff(o, NamedTuple{(:event_config,)}))
    config_for(o) = get(o, :event_config, event_config)

    shared_without = nothing
    if events_only
        println("\n######## Without events (shared by every value) ########")
        shared_without = run_scenario(config_for(overrides[1]); settings_for(overrides[1])...,
            run_validation = false, events_on = false,
            results_name = joinpath(name, "without_events"))
    end

    rows = []
    for (v, o) in zip(sweep_values, overrides)
        value_name = joinpath(name, "value_$v")
        println("\n######## $sweep_label = $v: with events ########")
        with = run_scenario(config_for(o); settings_for(o)..., run_validation,
                            results_name = joinpath(value_name, "with_events"))
        without = if events_only
            shared_without
        else
            println("\n######## $sweep_label = $v: without events ########")
            run_scenario(config_for(o); settings_for(o)..., run_validation = false, events_on = false,
                         results_name = joinpath(value_name, "without_events"))
        end
        df = compare_runs(with, without, normpath(joinpath(@__DIR__, "..", "Results", value_name)))
        push!(rows, summarise_comparison(v, df))
    end

    summary = DataFrame(identity.(rows))   # identity. gives the rows one concrete type
    write_sweep_summary(summary, sweep_label, events_only, gems_seed_used, folder)
    return (; summary, folder)
end


println("End Helpers")