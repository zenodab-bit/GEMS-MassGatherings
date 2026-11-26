# ===========================================================================
# events.jl
#
# Turns an EventConfig into concrete events: a date and size for every draw
# of every category, and for staged categories every venue of every stage.
# ===========================================================================

## === Sample Events ===

# For each Category: n_draws distinct dates within date_range, and one Event
# per section on each. Draws are numbered by date (draw 1 is the earliest).
# Dates are distinct per category only; categories sharing a date are
# handled by assign_events!.
function sample_events(event_config::EventConfig, rng)
    events = Vector{Event}()

    for category in event_config.categories
        if category isa StagedCategory
            append!(events, sample_staged_events(category, rng))
            continue
        end
        used_dates = Set{Int}()

        # random numbers per draw: the date, then each section's size;
        # then the draws are numbered by date
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

                # readable name from the parts that apply; the id if none do
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
                # short label: the section's name, else the category's
                base = !isempty(section.name) ? section.name :
                       isempty(category.name) ? id :
                       length(category.sections) > 1 ? "$(category.name)_sec$(section.id)" : category.name

                push!(events, Event(
                    id=id,
                    name=name,
                    label=with_draw(base, draw, category.n_draws),
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

    return unique_labels(events)
end

# A short label with the draw number when the category has several draws.
with_draw(base, draw, n_draws) = n_draws > 1 ? "$(base)_$(draw)" : base

# Labels shared by several events (e.g. two categories with a section
# "seated") are replaced by those events' full names, which are unique.
function unique_labels(events)
    labels = [e.label for e in events]
    clashes = Set(l for l in labels if count(==(l), labels) > 1)
    return [e.label in clashes ? reconstruct(e, label = e.name) : e for e in events]
end

# The same for a StagedCategory: one Event per venue of every stage, all on
# the draw's date. Random numbers per draw: the date, then each venue's
# size, stage by stage. section_id is the venue's number, stage_id the stage's.
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
            # short label: the venue's name, else the stage's, else the full name
            venue_part = length(stage.venues) > 1 ? "_venue$(v)" : ""
            base = !isempty(venue.name) ? venue.name :
                   !isempty(stage.name) ? stage.name * venue_part :
                   (isempty(category.name) ? "category$(category.id)" : category.name) * "_stage$(s)" * venue_part
            push!(events, Event(
                id="$(category.id)_$(draw)_$(s)_$(v)",
                name=join(parts, "_"),
                label=with_draw(base, draw, category.n_draws),
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

