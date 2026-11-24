# ===========================================================================
# config.jl
#
# The configuration structs a user builds in config/user_config.jl, the
# Event struct made from them, and validate_config.
# ===========================================================================

## === Configuration structs ===

# A group within a Category draw (e.g. "seated"), with its own attendees.
@with_kw struct EventSection
    id::Int
    name::String = ""
    n_range::Tuple{Int,Int}      # size, drawn uniformly from this range
    mean_contacts::Float64       # mean contacts per attendee within the section
    std_contacts::Float64 = 0.0  # Poisson instead of negative binomial when std^2 <= mean_contacts
end

# Supertype of Category and StagedCategory, so an EventConfig can hold both.
abstract type AbstractCategory end

# A kind of gathering (e.g. "festival") recurring n_draws times, each on
# its own date with its own attendees per section.
@with_kw struct Category <: AbstractCategory
    id::Int                                     # unique across categories
    name::String = ""                           # part of the event names, which must be unique
    date_range::Tuple{Int,Int}                  # (first_day, last_day) for the draws
    sections::Vector{EventSection}
    n_draws::Int                                # also the most draws one person can attend
    min_age::Int = 0                            # eligibility
    max_age::Int = 999
    age_weights::Vector{Float64} = Float64[]    # overrides EventConfig.age_dist if non-empty
    sex_weights::Vector{Float64} = Float64[]    # overrides EventConfig.sex_dist if non-empty
    core::Float64 = 0.0                         # share (0-1) of the smallest draw that attends every draw
    loyalty::Float64 = 0.0                      # share (0-1) of non-core spots kept for past non-core attendees
    cross_section_mean_contacts::Float64 = 0.0  # contacts with other sections of the draw; 0 keeps them apart
    cross_section_std_contacts::Float64 = 0.0
    # where core and loyal members go with several sections:
    #   :event   - a random section at every draw, in proportion to size
    #   :section - core keep a fixed home section; loyal return to their last one
    section_mode::Symbol = :event
end

# A place within a Stage (a restaurant, a stadium); no contacts with other venues.
@with_kw struct Venue
    name::String = ""
    n_range::Tuple{Int,Int}      # size, drawn anew at every draw
    mean_contacts::Float64       # mean contacts per attendee within the venue
    std_contacts::Float64 = 0.0  # Poisson instead of negative binomial when std^2 <= mean_contacts
end

# A stage of a StagedCategory: attendees drawn together, spread over its venues by size.
@with_kw struct Stage
    name::String = ""
    venues::Vector{Venue}
    overlap::Float64 = 0.0                      # share (0-1) taken from the previous stage's attendees; 0 for the first
    min_age::Int = 0                            # eligibility
    max_age::Int = 999
    age_weights::Vector{Float64} = Float64[]    # overrides EventConfig.age_dist if non-empty
    sex_weights::Vector{Float64} = Float64[]    # overrides EventConfig.sex_dist if non-empty
end

# Stages on the same day (e.g. restaurants and a match); people can attend
# several stages of a draw. No core or loyalty.
@with_kw struct StagedCategory <: AbstractCategory
    id::Int                                     # unique across categories
    name::String = ""                           # part of the event names, which must be unique
    date_range::Tuple{Int,Int}                  # (first_day, last_day) for the draws
    stages::Vector{Stage}                       # in order: overlap is with the previous one
    n_draws::Int
end

# The whole event setup.
@with_kw struct EventConfig
    categories::Vector{AbstractCategory}    # Category and/or StagedCategory
    age_boundaries::Vector{Int}              # increasing, e.g. [45, 65]
    age_dist::Vector{Float64}                # target age distribution; should sum to 1.0
    sex_dist::Vector{Vector{Float64}}        # target sex distribution per age group
end

# One concrete gathering on one date: a section of a Category draw, or a
# venue of a StagedCategory stage. Made by sample_events (events.jl).
@with_kw struct Event
    id::String                # "$(category)_$(draw)_$(section)", staged "$(category)_$(draw)_$(stage)_$(venue)"
    name::String               # readable; unique (validate_config checks)
    category_id::Int
    draw_id::Int                # 1..n_draws, in date order
    section_id::Int             # the section; for a staged category, the venue's number in its stage
    date::Int
    n::Int                      # target attendee count
    mean_contacts::Float64
    std_contacts::Float64
    cross_section_mean_contacts::Float64 = 0.0
    cross_section_std_contacts::Float64 = 0.0
    stage_id::Int = 0           # the stage (1, 2, ...); 0 for a normal Category
end


## === Validation ===

# Stops with an error on a configuration mistake (n_draws = 0 only warns).
# Called by config/user_config.jl, and by run_scenario before every run.
function validate_config(event_config::EventConfig)
    age_boundaries = event_config.age_boundaries
    isempty(age_boundaries) &&
        error("age_boundaries is empty; it needs at least one age, e.g. [45, 65].")
    issorted(age_boundaries) && allunique(age_boundaries) ||
        error("age_boundaries must be increasing without repeats, got $(age_boundaries).")
    n_groups = length(age_boundaries) + 1

    # event ids are built from category ids, so duplicates would merge events
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
        category.date_range[1] >= 0 ||
            error("Category $(category.name): date_range $(category.date_range) starts before day 0, " *
                  "the first simulated day.")

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
            section.n_range[1] >= 0 ||
                error("Category $(category.name), section $(section.id): n_range $(section.n_range) can't be negative.")
        end

        0.0 <= category.core <= 1.0 ||
            error("Category $(category.name): core ($(category.core)) must be between 0 and 1.")
        0.0 <= category.loyalty <= 1.0 ||
            error("Category $(category.name): loyalty ($(category.loyalty)) must be between 0 and 1.")

        category.section_mode in (:event, :section) ||
            error("Category $(category.name): section_mode must be :event or :section, got :$(category.section_mode).")
    end

    # per-event plots are named after events; names don't depend on the rng
    names = [e.name for e in sample_events(event_config, Xoshiro(0))]
    duplicates = unique(n for n in names if count(==(n), names) > 1)
    isempty(duplicates) ||
        error("Several events would be named $(join(repr.(duplicates), ", ")), so their plots would " *
              "overwrite each other: give the categories (or their sections, stages or venues) different names.")

    println("Config validated OK.")
end

# The same for a StagedCategory's stages and venues.
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
            venue.n_range[1] >= 0 ||
                error("$where, venue $v ($(venue.name)): n_range $(venue.n_range) can't be negative.")
        end
    end
end

