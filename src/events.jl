# ===========================================================================
# events.jl
#
# Turns an EventConfig into concrete events: dates and sizes for every draw
# of every category (sample_events), and for staged categories every venue
# of every stage (sample_staged_events).
# ===========================================================================

## === Sample Events ===

# Turns an EventConfig into concrete Event objects: for each Category,
# draws n_draws distinct dates within its date_range, and for each date
# creates one Event per EventSection (sections of one draw share that date).
# Draws are numbered by date: draw 1 is the earliest.
#
# Date uniqueness is tracked per category, not globally, so two different
# categories can coincidentally land on the same date;
# the same-day conflict check in assign_events! handles it.
function sample_events(event_config::EventConfig, rng)
    events = Vector{Event}()

    for category in event_config.categories
        if category isa StagedCategory
            append!(events, sample_staged_events(category, rng))
            continue
        end
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

# The same for a StagedCategory: n_draws distinct dates within its
# date_range, each draw with one Event per venue of every stage, all on the
# draw's date. Per draw, the random numbers are taken in this order: the
# date, then each venue's size, stage by stage. Event ids are
# "$(category)_$(draw)_$(stage)_$(venue)"; section_id holds the venue's
# number within its stage, stage_id the stage's.
function sample_staged_events(category::StagedCategory, rng)
    used_dates = Set{Int}()
    draws = Tuple{Int,Vector{Vector{Int}}}[]   # (date, size per venue, per stage)
    for _ in 1:category.n_draws
        available_dates = setdiff(category.date_range[1]:category.date_range[2], used_dates)
        isempty(available_dates) &&
            error("Staged category $(category.name): date_range $(category.date_range) has only " *
                  "$(length(category.date_range[1]:category.date_range[2])) day(s) available, " *
                  "not enough unique dates for $(category.n_draws) draws.")
        date = rand(rng, available_dates)
        push!(used_dates, date)
        push!(draws, (date, [[rand(rng, v.n_range[1]:v.n_range[2]) for v in stage.venues]
                             for stage in category.stages]))
    end
    sort!(draws, by=first)   # draw 1 is the earliest

    events = Vector{Event}()
    for (draw, (date, sizes)) in enumerate(draws)
        for (s, stage) in enumerate(category.stages), (v, venue) in enumerate(stage.venues)
            parts = String[]
            isempty(category.name) || push!(parts, category.name)
            category.n_draws > 1 && push!(parts, "draw$(draw)")
            push!(parts, isempty(stage.name) ? "stage$(s)" : stage.name)
            length(stage.venues) > 1 && push!(parts, isempty(venue.name) ? "venue$(v)" : venue.name)
            push!(events, Event(
                id="$(category.id)_$(draw)_$(s)_$(v)",
                name=join(parts, "_"),
                category_id=category.id,
                draw_id=draw,
                section_id=v,
                date=date,
                n=sizes[s][v],
                mean_contacts=venue.mean_contacts,
                std_contacts=venue.std_contacts,
                stage_id=s
            ))
        end
    end
    return events
end


