# ===========================================================================
# config.jl
#
# The event configuration: the structs a user builds in config/user_config.jl
# (Section, Category, Venue, Stage, StagedCategory, EventConfig), the Event
# struct the code creates from them, and the checks of a configuration
# (validate_config).
# ===========================================================================

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
# Produced by sample_events() (events.jl) — not user-constructed.
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
# construction (see end of config/user_config.jl).
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


println("End config")
