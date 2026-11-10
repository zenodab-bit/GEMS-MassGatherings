# ===========================================================================
# 3_Population.jl
#
# Turns an EventConfig into concrete events and a simulation-ready
# population, and assigns people to those events. Three pieces, used
# together by 2_Interface.jl in this order:
#   1. sample_events(event_config, rng)       - returns Vector{Event}
#   2. prepare_population(event_config, rng, people, transmission_dist) - returns people DataFrame + labels
#   3. assign_events!(people, events, ...)    - returns attendees Dict and core_groups,
#                                                mutates `people` in place

# ===========================================================================


## === Sample Events ===

# Turns an EventConfig into concrete Event objects: for each Category,
# draws n_draws distinct dates within its date_range, and for each date
# creates one Event per Section (sections of one draw share that date).
# Draws are numbered by date: draw 1 is the earliest.
#
# Date uniqueness is tracked per category, not globally, so two different
# categories can coincidentally land on the same date;
# the same-day conflict check in assign_events! handles it.
function sample_events(event_config::EventConfig, rng)
    events = Vector{Event}()

    for category in event_config.categories
        used_dates = Set{Int}()  # track dates used by this category

        # sample every draw first (same RNG order as before: date, then each
        # section's size), then number draws by date so draw 1 is the earliest
        draws = Tuple{Int,Vector{Int}}[]   # (date, size per section)
        for _ in 1:category.n_draws
            available_dates = setdiff(category.date_range[1]:category.date_range[2], used_dates)

            if isempty(available_dates)
                error("Category $(category.name): date_range $(category.date_range) has only " *
                      "$(length(category.date_range[1]:category.date_range[2])) day(s) available, " *
                      "not enough unique dates for $(category.n_draws) draws.")
            end

            date = rand(rng, available_dates)
            push!(used_dates, date)
            push!(draws, (date, [rand(rng, s.n_range[1]:s.n_range[2]) for s in category.sections]))
        end
        sort!(draws, by=first)

        for (draw, (date, sizes)) in enumerate(draws)
            for (section, n) in zip(category.sections, sizes)
                id = "$(category.id)_$(draw)_$(section.id)"

                # Build a readable name from whatever parts apply; falls back to id.
                parts = String[]
                if !isempty(category.name)
                    push!(parts, category.name)
                end
                if category.n_draws > 1
                    push!(parts, "draw$(draw)")
                end
                if length(category.sections) > 1
                    if !isempty(section.name)
                        push!(parts, section.name)
                    else
                        push!(parts, "sec$(section.id)")
                    end
                end
                name = isempty(parts) ? id : join(parts, "_")

                push!(events, Event(
                    id=id,
                    name=name,
                    category_id=category.id,
                    draw_id=draw,
                    section_id=section.id,
                    date=date,
                    n=n,
                    mean_contacts=section.mean_contacts,
                    std_contacts=section.std_contacts,
                    cross_section_mean_contacts=category.cross_section_mean_contacts,
                    cross_section_std_contacts=category.cross_section_std_contacts
                ))
            end
        end
    end

    return events
end


## === Prepare Population ===

# Takes the base population table (chosen in 2_Interface.jl, or the
# synthetic one in test.jl) and adds every column the rest of the project
# needs: age-group labels, empty per-person event-attendance records, and
# each person's transmission_prob. The table is changed in place.
#
# `people` needs at least id, age and sex (1 or 2) columns, plus the GEMS
# setting ids (household, and optionally office, schoolclass, municipality).
#
# transmission_prob is drawn from `transmission_dist` (set in
# 1_UserConfig.jl) and kept within [0, 1], since it's used as a probability.
#
# The attendance columns added here (category_ids, draw_ids, section_ids,
# mean_event_contacts, std_event_contacts, event_dates) are also listed in
# 2_Interface.jl's ind_extension, which is how GEMS carries them onto its
# own Individual objects for 4_Contacts.jl and 5_Transmission.jl to read.
function prepare_population(event_config::EventConfig, rng, people::DataFrame, transmission_dist)

    age_boundaries = event_config.age_boundaries
    age_groups = [age_group_label_from_idx(i, age_boundaries) for i in 1:(length(age_boundaries)+1)]
    sex_levels = [1, 2]

    people.age_group = age_group_label.(people.age, Ref(age_boundaries))
    people.age_group = categorical(people.age_group; ordered=true, levels=age_groups)

    people.category_ids = [Int32[] for _ in 1:nrow(people)]
    people.draw_ids = [Int32[] for _ in 1:nrow(people)]
    people.section_ids = [Int32[] for _ in 1:nrow(people)]
    people.mean_event_contacts = [Float64[] for _ in 1:nrow(people)]
    people.event_dates = [Int32[] for _ in 1:nrow(people)]
    people.std_event_contacts = [Float64[] for _ in 1:nrow(people)]
    people.cross_section_mean_contacts = [Float64[] for _ in 1:nrow(people)]
    people.cross_section_std_contacts = [Float64[] for _ in 1:nrow(people)]

    # clamped, since e.g. a Gamma can go above 1 and a Normal below 0
    people.transmission_prob = Float64[clamp(rand(rng, transmission_dist), 0.0, 1.0) for _ in 1:nrow(people)]

    return people, age_groups, sex_levels
end


## === Compute Weights ===

# Computes per-candidate sampling weights for demographic-targeted event
# selection, used by assign_events! for both core-group and random-pool draws.
#
# Weights are per-capita: each candidate's weight is the target share
# (age_weight times sex_weight) divided by how many candidates share their
# (age_group, sex) cell. This keeps a cell's total pull equal to its
# intended target share regardless of how many candidates happen to fall
# into it — without dividing, a larger candidate pool would pull more than
# its intended share just by having more people in it.
#
# Cells with zero candidates get weight 0; if their target share is
# nonzero, sampling redistributes it proportionally across the other cells.
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

        # category age_weights/sex_weights override EventConfig.age_dist/sex_dist if given
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

# Splits k items across slots in proportion to `sizes`, using largest
# remainders so the counts always add up to exactly k.
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

# Demographically-weighted draw of up to k people from pool, without
# replacement. Returns an empty vector, without touching rng, when there's
# nothing to draw.
function weighted_draw(rng, people, pool, k, category, event_config)
    k = min(k, length(pool))
    k <= 0 && return Int[]
    weights = compute_weights(people, pool, category, event_config)
    return sample(rng, pool, Weights(weights), k, replace=false)
end

# Hands out `members` to the events in `evs`, in consecutive blocks of
# `counts`, appending to dest[event.id].
function split_into!(dest::Dict{String,Vector{Int}}, evs, members, counts)
    pos = 0
    for (e, c) in zip(evs, counts)
        append!(dest[e.id], members[pos+1:pos+c])
        pos += c
    end
    return dest
end

# Assigns people to every event, mutating `people`'s attendance columns and
# returning attendees::Dict{String,Vector{Int}} (event id to attendee
# row-indices) and core_groups::Dict{Int,Vector{Int}} (category id to core member row-indices)
#
# Events don't need to be pre-sorted: dates are processed in calendar order
# internally, which loyalty tracking depends on.
#
# Phase 1 builds each category's core group once (a person can be core for
# at most one category). Phase 2 goes date by date and fills all events of a
# date together, in three rounds: every event's core, then every event's
# loyalty spots, then every event's random fill. On a shared date, priority
# goes by strength of claim (core > loyalty > random), not by category id.
# A person still attends at most one event per day.
#
# With several sections, section_mode decides where core and loyal members
# go: :event spreads them over sections at every draw, in proportion to
# size; :section gives core a fixed home section and sends loyal attendees
# back to their last section.

function assign_events!(people::DataFrame, events::Vector{Event}, event_config::EventConfig, rng)
    loyal_pools = Dict{Int,Set{Int}}()
    core_groups = Dict{Int,Vector{Int}}()
    all_core_members = Set{Int}()
    attendees = Dict{String,Vector{Int}}()
    base_candidates_by_category = Dict{Int,Vector{Int}}()

    categories_by_id = Dict(c.id => c for c in event_config.categories)
    home_section = Dict{Int,Int}()              # core member => section_id (:section mode only)
    last_section = Dict{Int,Dict{Int,Int}}()    # category_id => (person => section last attended)

    # --- PHASE 1: build each category's core group once ---
    for category in event_config.categories
        cat_events = filter(e -> e.category_id == category.id, events)
        isempty(cat_events) && continue

        # core is a fraction of the whole draw (all sections together).
        # :event uses the smallest draw; :section uses each section's smallest
        # size, so every home section's share always fits.
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

        # demographically weighted, excluding anyone already core elsewhere
        if n_core_total > 0
            remaining_for_core = filter(idx -> idx ∉ all_core_members, base_candidates)
            core_weights = compute_weights(people, remaining_for_core, category, event_config)
            n_to_sample = min(n_core_total, length(remaining_for_core))
            append!(core_groups[category.id],
                sample(rng, remaining_for_core, Weights(core_weights), n_to_sample, replace=false))
        end

        union!(all_core_members, core_groups[category.id])

        # :section mode: fixed home section per core member, in proportion to
        # each section's smallest size. Each section's share is then at most
        # its smallest size, so a home section can never be overfilled.
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

    # --- PHASE 2: fill events date by date, in three rounds per date ---
    for date in sort!(unique(e.date for e in events))
        day_events = sort(filter(e -> e.date == date, events),
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
                # first: back to the section they last attended
                taken = Set{Int}()
                for (k, e) in enumerate(secs)
                    own = [idx for idx in eligible if get(last_sec, idx, 0) == e.section_id]
                    picked = weighted_draw(rng, people, own, quotas[k], category, event_config)
                    append!(loyalty_selected[e.id], picked)
                    union!(taken, picked)
                end
                # then: open spots go to loyal people whose own section was full
                rest = [idx for idx in eligible if idx ∉ taken]
                for (k, e) in enumerate(secs)
                    gap = quotas[k] - length(loyalty_selected[e.id])
                    picked = weighted_draw(rng, people, rest, gap, category, event_config)
                    append!(loyalty_selected[e.id], picked)
                    picked_set = Set(picked)
                    filter!(idx -> idx ∉ picked_set, rest)
                end
            else
                # one draw for the whole event, shuffled, then split by quota
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

        # Round 3: random fill for whatever is left
        for e in day_events
            category = categories_by_id[e.category_id]
            pool = filter(idx -> idx ∉ booked_today, base_candidates_by_category[e.category_id])
            picked = weighted_draw(rng, people, pool, e.n - length(selected[e.id]), category, event_config)
            random_selected[e.id] = picked
            append!(selected[e.id], picked)
            union!(booked_today, picked)
        end

        # Record attendance; loyalty/random attendees become loyalty-eligible
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
            end
        end
    end
    return attendees, core_groups
end


println("End Population")