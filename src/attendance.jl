# ===========================================================================
# attendance.jl
#
# Prepares the population (event columns, transmission_prob) and assigns
# people to events: core, loyalty and random fill for categories, overlap
# and fill stage by stage for staged categories.
# ===========================================================================

## === Prepare Population ===

# The per-person event columns (name => element type): created empty by
# prepare_population, filled by assign_events!, copied by GEMS onto its
# individuals for contacts.jl and transmission.jl.
const EVENT_COLUMN_TYPES = [
    :category_ids => Int32,
    :draw_ids => Int32,
    :section_ids => Int32,
    :mean_event_contacts => Float64,
    :event_dates => Int32,
    :std_event_contacts => Float64,
    :cross_section_mean_contacts => Float64,
    :cross_section_std_contacts => Float64,
    :stage_ids => Int32,   # stage of a staged category; 0 for a normal one
]

# All columns GEMS needs on its individuals (ind_extension in run_scenario).
const GEMS_COLUMNS = [first.(EVENT_COLUMN_TYPES); :transmission_prob]

# Adds to `people` (in place) age groups, the empty event columns and each
# person's transmission_prob.
function prepare_population(event_config::EventConfig, rng, people::DataFrame, transmission_dist)

    age_boundaries = event_config.age_boundaries
    age_groups = [age_group_label_from_idx(i, age_boundaries) for i in 1:(length(age_boundaries)+1)]
    sex_levels = [1, 2]

    people.age_group = age_group_label.(people.age, Ref(age_boundaries))
    people.age_group = categorical(people.age_group; ordered=true, levels=age_groups)

    for (col, T) in EVENT_COLUMN_TYPES
        people[!, col] = [T[] for _ in 1:nrow(people)]
    end

    # clamped: e.g. a Gamma can exceed 1, a Normal go below 0
    people.transmission_prob = Float64[clamp(rand(rng, transmission_dist), 0.0, 1.0) for _ in 1:nrow(people)]

    return people, age_groups, sex_levels
end


## === Compute Weights ===

# Per candidate: the target share of their (age group, sex) cell divided by
# the cell's candidates, so each cell pulls its target share.
function compute_weights(people, candidates, category, event_config)
    age_boundaries = event_config.age_boundaries
    n_groups = length(age_boundaries) + 1

    ages = people.age
    sexes = people.sex

    gidx = Vector{Int}(undef, length(candidates))
    for (k, idx) in enumerate(candidates)
        gidx[k] = age_group_idx(ages[idx], age_boundaries)
    end

    counts = zeros(Int, n_groups, 2)
    for (k, idx) in enumerate(candidates)
        counts[gidx[k], sexes[idx]] += 1
    end

    weights = Vector{Float64}(undef, length(candidates))
    for (k, idx) in enumerate(candidates)
        g = gidx[k]
        sex = sexes[idx]

        # the category's (or stage's) weights override age_dist / sex_dist
        aw = isempty(category.age_weights) ?
             event_config.age_dist[g] :
             category.age_weights[g]

        sw = isempty(category.sex_weights) ?
             event_config.sex_dist[g][sex] :
             category.sex_weights[sex]

        c = counts[g, sex]
        weights[k] = c > 0 ? (aw * sw) / c : 0.0
    end

    return weights
end


## === Assign Events ===

# Splits k across slots in proportion to `sizes` (largest remainders, so
# the counts add up to exactly k).
function proportional_counts(k::Integer, sizes::AbstractVector{<:Real})
    counts = zeros(Int, length(sizes))
    (k <= 0 || isempty(sizes)) && return counts
    total = sum(sizes)
    if total <= 0
        counts[1] = k
        return counts
    end
    raw = k .* sizes ./ total
    counts .= floor.(Int, raw)
    order = sortperm(raw .- counts, rev=true)
    for i in 1:(k-sum(counts))
        counts[order[i]] += 1
    end
    return counts
end

# Weighted draw of up to k people from pool, without replacement. People
# with weight 0 are never drawn, so it can return fewer than k (callers
# warn about unfilled places). Doesn't touch rng when there's nothing to draw.
function weighted_draw(rng, people, pool, k, category, event_config)
    (k <= 0 || isempty(pool)) && return Int[]
    weights = compute_weights(people, pool, category, event_config)
    k = min(k, count(>(0), weights))
    k <= 0 && return Int[]
    return sample(rng, pool, Weights(weights), k, replace=false)
end

# Hands `members` out to the events in `evs` in consecutive blocks of
# `counts`, appending to dest[event.id].
function split_into!(dest::Dict{String,Vector{Int}}, evs, members, counts)
    pos = 0
    for (e, c) in zip(evs, counts)
        append!(dest[e.id], members[pos+1:pos+c])
        pos += c
    end
    return dest
end

# Assigns people to every event; returns attendees (event id => rows) and
# core_groups (category id => rows). Phase 1 picks each core group once.
# Phase 2 goes date by date in calendar order (the loyal pools depend on
# it), in three rounds over all of a date's events: cores, loyalty, random
# fill. So on a shared date the stronger claim wins, whatever the category
# id. Staged categories come last. Nobody attends two events on a day.
function assign_events!(people::DataFrame, events::Vector{Event}, event_config::EventConfig, rng)
    loyal_pools = Dict{Int,Set{Int}}()
    core_groups = Dict{Int,Vector{Int}}()
    all_core_members = Set{Int}()
    attendees = Dict{String,Vector{Int}}()
    base_candidates_by_category = Dict{Int,Vector{Int}}()

    categories_by_id = Dict(c.id => c for c in event_config.categories)
    home_section = Dict{Int,Int}()              # core member => section_id (:section mode only)
    last_section = Dict{Int,Dict{Int,Int}}()    # category_id => (person => section last attended)

    # --- PHASE 1: each category's core group ---
    for category in event_config.categories
        category isa StagedCategory && continue   # no core in staged categories
        cat_events = filter(e -> e.category_id == category.id, events)
        isempty(cat_events) && continue

        # core is a share of the smallest draw (:event), or of each
        # section's smallest size (:section), so it always fits
        section_ids = [s.id for s in category.sections]
        min_section_n = [minimum(e.n for e in cat_events if e.section_id == sid) for sid in section_ids]
        if category.section_mode == :section
            core_base = sum(min_section_n)
        else
            draw_totals = [sum(e.n for e in cat_events if e.draw_id == d)
                           for d in unique(e.draw_id for e in cat_events)]
            core_base = minimum(draw_totals)
        end
        n_core_total = round(Int, core_base * category.core)

        base_candidates = findall(
            (people.age .>= category.min_age) .&
            (people.age .<= category.max_age)
        )
        base_candidates_by_category[category.id] = base_candidates

        core_groups[category.id] = Int[]

        if n_core_total > 0
            remaining_for_core = filter(idx -> idx ∉ all_core_members, base_candidates)
            append!(core_groups[category.id],
                weighted_draw(rng, people, remaining_for_core, n_core_total, category, event_config))
            length(core_groups[category.id]) < n_core_total &&
                @warn "Category $(category.name): core group of $n_core_total wanted, but only " *
                      "$(length(core_groups[category.id])) eligible people (with a weight above 0, " *
                      "not core elsewhere)."
        end

        union!(all_core_members, core_groups[category.id])

        # :section mode: a fixed home section per core member, in proportion
        # to each section's smallest size, so no home section overfills
        if category.section_mode == :section && length(category.sections) > 1
            members = shuffle(rng, core_groups[category.id])
            counts = proportional_counts(length(members), min_section_n)
            pos = 0
            for (sid, c) in zip(section_ids, counts)
                for idx in members[pos+1:pos+c]
                    home_section[idx] = sid
                end
                pos += c
            end
        end
    end

    # --- PHASE 2: date by date, three rounds per date ---
    for date in sort!(unique(e.date for e in events))
        # normal categories' events (staged ones after them, below)
        day_events = sort(filter(e -> e.date == date && e.stage_id == 0, events),
            by=e -> (e.category_id, e.draw_id, e.section_id))
        # one draw per category per date: (category, its section events)
        day_draws = [(categories_by_id[cid], filter(e -> e.category_id == cid, day_events))
                     for cid in unique(e.category_id for e in day_events)]

        booked_today = Set{Int}()
        selected = Dict(e.id => Int[] for e in day_events)
        loyalty_selected = Dict(e.id => Int[] for e in day_events)
        random_selected = Dict(e.id => Int[] for e in day_events)

        # Round 1: core
        for (category, secs) in day_draws
            core = filter(idx -> idx ∉ booked_today, get(core_groups, category.id, Int[]))
            isempty(core) && continue
            if length(secs) == 1
                append!(selected[secs[1].id], core)
            elseif category.section_mode == :section
                sec_event = Dict(e.section_id => e for e in secs)
                for idx in core
                    push!(selected[sec_event[home_section[idx]].id], idx)
                end
            else  # :event
                members = shuffle(rng, core)
                counts = proportional_counts(length(members), [e.n for e in secs])
                split_into!(selected, secs, members, counts)
            end
            union!(booked_today, core)
        end

        # Round 2: loyalty (quota per section = loyalty × spots left after core)
        for (category, secs) in day_draws
            (category.loyalty > 0.0 && haskey(loyal_pools, category.id)) || continue
            quotas = Int[]
            for e in secs
                spots = max(0, e.n - length(selected[e.id]))
                push!(quotas, min(round(Int, category.loyalty * spots), spots))
            end
            eligible = [idx for idx in loyal_pools[category.id] if idx ∉ booked_today]

            if category.section_mode == :section && length(secs) > 1
                last_sec = get!(() -> Dict{Int,Int}(), last_section, category.id)
                # first back to their last section
                taken = Set{Int}()
                for (k, e) in enumerate(secs)
                    own = [idx for idx in eligible if get(last_sec, idx, 0) == e.section_id]
                    picked = weighted_draw(rng, people, own, quotas[k], category, event_config)
                    append!(loyalty_selected[e.id], picked)
                    union!(taken, picked)
                end
                # then open spots to those whose own section was full
                rest = [idx for idx in eligible if idx ∉ taken]
                for (k, e) in enumerate(secs)
                    gap = quotas[k] - length(loyalty_selected[e.id])
                    picked = weighted_draw(rng, people, rest, gap, category, event_config)
                    append!(loyalty_selected[e.id], picked)
                    picked_set = Set(picked)
                    filter!(idx -> idx ∉ picked_set, rest)
                end
            else
                # one draw for the whole event, shuffled, split by quota
                picked = weighted_draw(rng, people, eligible, sum(quotas), category, event_config)
                if length(secs) > 1
                    picked = shuffle(rng, picked)
                end
                split_into!(loyalty_selected, secs, picked, proportional_counts(length(picked), quotas))
            end

            for e in secs
                append!(selected[e.id], loyalty_selected[e.id])
                union!(booked_today, loyalty_selected[e.id])
            end
        end

        # Round 3: random fill
        for e in day_events
            category = categories_by_id[e.category_id]
            pool = filter(idx -> idx ∉ booked_today, base_candidates_by_category[e.category_id])
            picked = weighted_draw(rng, people, pool, e.n - length(selected[e.id]), category, event_config)
            random_selected[e.id] = picked
            append!(selected[e.id], picked)
            union!(booked_today, picked)
            length(selected[e.id]) < e.n &&
                @warn "Event $(e.id) on day $(e.date): only $(length(selected[e.id])) of $(e.n) places filled; " *
                      "not enough eligible people left that day."
        end

        # record attendance; loyalty and random attendees join the loyal pool
        for e in day_events
            pool = get!(() -> Set{Int}(), loyal_pools, e.category_id)
            last_sec = get!(() -> Dict{Int,Int}(), last_section, e.category_id)
            for idx in Iterators.flatten((loyalty_selected[e.id], random_selected[e.id]))
                push!(pool, idx)
                last_sec[idx] = e.section_id
            end

            attendees[e.id] = selected[e.id]
            for idx in selected[e.id]
                push!(people.category_ids[idx], Int32(e.category_id))
                push!(people.draw_ids[idx], Int32(e.draw_id))
                push!(people.section_ids[idx], Int32(e.section_id))
                push!(people.mean_event_contacts[idx], e.mean_contacts)
                push!(people.std_event_contacts[idx], e.std_contacts)
                push!(people.event_dates[idx], Int32(e.date))
                push!(people.cross_section_mean_contacts[idx], e.cross_section_mean_contacts)
                push!(people.cross_section_std_contacts[idx], e.cross_section_std_contacts)
                push!(people.stage_ids[idx], Int32(0))
            end
        end

        # staged categories, after every other category of the day
        staged_today = sort(filter(e -> e.date == date && e.stage_id > 0, events),
            by=e -> (e.category_id, e.stage_id, e.section_id))
        for cid in unique(e.category_id for e in staged_today)
            fill_staged_draw!(attendees, people, categories_by_id[cid],
                filter(e -> e.category_id == cid, staged_today), booked_today, event_config, rng)
        end
    end
    return attendees, core_groups
end

# Fills one staged draw, stage by stage: the overlap from the previous
# stage, the rest from people not yet booked that day; then spread over the
# venues by size. Warns about places it can't fill.
function fill_staged_draw!(attendees, people, category::StagedCategory, draw_events,
                           booked_today::Set{Int}, event_config, rng)
    in_draw = Set{Int}()      # attendees of this draw's stages so far
    previous = Int[]          # the previous stage's attendees
    for (s, stage) in enumerate(category.stages)
        venues = sort(filter(e -> e.stage_id == s, draw_events), by=e -> e.section_id)
        total = sum(e.n for e in venues)
        eligible(idx) = stage.min_age <= people.age[idx] <= stage.max_age

        n_overlap = s == 1 ? 0 : round(Int, stage.overlap * total)
        from_previous = weighted_draw(rng, people, filter(eligible, previous), n_overlap, stage, event_config)
        length(from_previous) < n_overlap &&
            @warn "Staged category $(category.name), day $(venues[1].date), stage $s: overlap of $n_overlap " *
                  "wanted, but only $(length(from_previous)) eligible people at the previous stage."

        pool = filter(idx -> eligible(idx) && idx ∉ booked_today && idx ∉ in_draw, 1:nrow(people))
        fresh = weighted_draw(rng, people, pool, total - length(from_previous), stage, event_config)

        members = shuffle(rng, vcat(from_previous, fresh))
        counts = proportional_counts(length(members), [e.n for e in venues])
        pos = 0
        for (e, c) in zip(venues, counts)
            attendees[e.id] = members[pos+1:pos+c]
            pos += c
            c < e.n && @warn "Event $(e.id) on day $(e.date): only $c of $(e.n) places filled; " *
                              "not enough eligible people left that day."
            for idx in attendees[e.id]
                push!(people.category_ids[idx], Int32(e.category_id))
                push!(people.draw_ids[idx], Int32(e.draw_id))
                push!(people.section_ids[idx], Int32(e.section_id))
                push!(people.mean_event_contacts[idx], e.mean_contacts)
                push!(people.std_event_contacts[idx], e.std_contacts)
                push!(people.event_dates[idx], Int32(e.date))
                push!(people.cross_section_mean_contacts[idx], 0.0)
                push!(people.cross_section_std_contacts[idx], 0.0)
                push!(people.stage_ids[idx], Int32(s))
            end
        end
        union!(in_draw, members)
        previous = members
    end
    union!(booked_today, in_draw)
    return attendees
end

