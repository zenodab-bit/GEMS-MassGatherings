# ===========================================================================
# 7_Validation.jl
#
# Auto-run at the end of 2_Interface.jl, after 6_Analysis.jl, if
# run_validation is true. Five independent checks against the finished
# simulation: did assignment hit its targets, did core attendance per
# category work as configured, do compartment totals add up, do
# demographics match the config, and does the epidemic state roughly
# match what the population-wide infection log implies.
# ===========================================================================

# Population-wide SEIR compartment counts as of `tick`, computed directly
# from the infection log (not from tick_cases, which turned out to report
# same-day incidence rather than cumulative prevalence). Used by
# validate_epidemic_state as the population-wide baseline.
function population_compartment_counts(inf_log::DataFrame, tick::Int, pop_size::Int)
    infectious = 0
    exposed = 0
    recovered = 0
    dead = 0
    same_day_other = 0   # infected on `tick` outside GlobalSetting, before events run

    for row in eachrow(inf_log)
        if row.tick == tick
            row.setting_type != 'g' && (same_day_other += 1)
            continue
        end
        row.tick > tick && continue
        if row.infectiousness_onset <= tick &&
           (row.recovery > tick || row.recovery == -1) &&
           (row.death > tick || row.death == -1)
            infectious += 1
        elseif row.recovery != -1 && row.recovery <= tick
            recovered += 1
        elseif row.death != -1 && row.death <= tick
            dead += 1
        else
            exposed += 1
        end
    end

    susceptible = pop_size - infectious - exposed - recovered - dead - same_day_other
    return (susceptible=susceptible, infectious=infectious, exposed=exposed,
        recovered=recovered, dead=dead, same_day_other=same_day_other)
end

# Prints the overall attendance distribution (events attended per person,
# all categories together), checks that nobody is booked twice on the same
# day, and confirms each event hit its target attendee count.
# Core/non-core breakdowns are per category, in validate_category_attendance.
function validate_assignment(people::DataFrame, events::Vector{Event}, event_config::EventConfig, attendees::Dict{String,Vector{Int}})
    println("\n=== Assignment Validation ===")

    attendance_per_person = Dict{Int,Int}()
    for idxs in values(attendees)
        for idx in idxs
            attendance_per_person[idx] = get(attendance_per_person, idx, 0) + 1
        end
    end
    total_assigned = length(attendance_per_person)
    multi_assigned = count(>(1), values(attendance_per_person))

    println("\nAttendance distribution (all categories together):")
    println("  0 events: $(nrow(people) - total_assigned) people (unassigned)")
    attendance_counts = countmap(values(attendance_per_person))
    for n in sort(collect(keys(attendance_counts)))
        println("  $n event(s): $(attendance_counts[n]) people")
    end

    println("\nPeople attending at least one event: $total_assigned")
    println("People attending multiple events:    $multi_assigned")

    # nobody may attend two events on the same day
    n_double = count(d -> !allunique(d), people.event_dates)
    println(n_double == 0 ? "✓ no same-day double bookings" :
                            "✗ $n_double people booked twice on the same day")

    for event in events
        assigned = length(attendees[event.id])
        status = assigned == event.n ? "✓" : "✗ expected $(event.n)"
        println("Event $(event.id) (date=$(event.date), n=$(event.n)): $assigned assigned $status")
    end

    println("Unassigned: $(nrow(people) - total_assigned) / $(nrow(people))")
end

# Per category: how many people attended 1..n_draws draws of that category
# (core vs non-core), whether every core member attended every draw, and,
# with several sections, how many core members went to each section.
# Uses the real core groups returned by assign_events!, not a recomputed size.
function validate_category_attendance(events::Vector{Event}, event_config::EventConfig,
    attendees::Dict{String,Vector{Int}}, core_groups::Dict{Int,Vector{Int}})
    println("\n=== Per-category Attendance ===")

    for category in event_config.categories
        cat_events = filter(e -> e.category_id == category.id, events)
        isempty(cat_events) && continue

        # draws of this category attended, per person
        draws_per_person = Dict{Int,Int}()
        for e in cat_events, idx in attendees[e.id]
            draws_per_person[idx] = get(draws_per_person, idx, 0) + 1
        end

        core = Set(get(core_groups, category.id, Int[]))
        println("\n  $(category.name) (id=$(category.id), n_draws=$(category.n_draws), core=$(length(core))):")

        for k in sort(unique(values(draws_per_person)))
            people_k = [idx for (idx, c) in draws_per_person if c == k]
            n_core = count(in(core), people_k)
            println("    $k draw(s): $(length(people_k)) people ($n_core core + $(length(people_k) - n_core) non-core)")
        end

        # core guarantee: every core member attends every draw
        if isempty(core)
            println("    (no core)")
        else
            n_missed = count(idx -> get(draws_per_person, idx, 0) < category.n_draws, core)
            println(n_missed == 0 ?
                    "    ✓ all $(length(core)) core members attended every draw" :
                    "    ✗ $n_missed core member(s) missed at least one draw")
        end

        # core per section (multi-section categories only)
        if !isempty(core) && length(category.sections) > 1
            println("    Core per section:")
            for e in sort(cat_events, by=e -> (e.draw_id, e.section_id))
                println("      $(e.name): $(count(in(core), attendees[e.id])) core of $(length(attendees[e.id]))")
            end
        end
    end
end

# Checks, per event: infected_at_event never exceeds the event size, and
# every compartment sums back to the event size (within rounding).
function validate_results(aggregated::Dict, events::Vector{Event})
    println("\n=== Results Validation ===")

    for event in events
        metrics = aggregated[event.id]
        println("\nEvent $(event.id) — Day $(event.date) — n=$(event.n)")

        # infected at event should not exceed section size
        if metrics[:infected_at_event].max > event.n
            println("  ✗ infected_at_event max ($(metrics[:infected_at_event].max)) exceeds section size ($(event.n))")
        else
            println("  ✓ infected_at_event within bounds")
        end

        # susceptible + infected + exposed + recovered + dead should equal section size
        total = metrics[:susceptible].mean + metrics[:infectious].mean +
                metrics[:exposed].mean + metrics[:recovered].mean +
                metrics[:dead].mean + metrics[:same_day_other].mean
        if abs(total - event.n) > 1.0  # allow small rounding error
            println("  ✗ compartment totals ($(fmt(total))) don't match section size ($(event.n))")
        else
            println("  ✓ compartment totals match section size")
        end
    end
end


# Compares each event's actual attendee age/sex distribution against the
# category's (or, if unset, the population-wide) target.
function validate_demographics(people::DataFrame, events::Vector{Event}, event_config::EventConfig, attendees::Dict{String,Vector{Int}})
    println("\n=== Demographic Validation ===")

    for category in event_config.categories
        cat_events = filter(e -> e.category_id == category.id, events)
        isempty(cat_events) && continue

        # expected share per (age group, sex) cell, built the same way
        # compute_weights builds sampling weights: age weight x sex weight
        n_groups = length(event_config.age_boundaries) + 1
        mass = [(isempty(category.age_weights) ? event_config.age_dist[g] : category.age_weights[g]) *
                (isempty(category.sex_weights) ? event_config.sex_dist[g][s] : category.sex_weights[s])
                for g in 1:n_groups, s in 1:2]
        mass ./= sum(mass)
        expected_age = vec(sum(mass, dims=2))
        expected_sex = vec(sum(mass, dims=1))

        println("\nCategory $(category.id) — loyalty=$(category.loyalty):")

        for event in cat_events
            idxs = attendees[event.id]
            n = length(idxs)
            ages = people.age[idxs]
            sexes = people.sex[idxs]

            println("\n  Event $(event.id) (day=$(event.date), n=$(event.n)):")

            println("    Age distribution:")
            for i in 1:length(event_config.age_boundaries)
                boundary = event_config.age_boundaries[i]
                expected = expected_age[i]
                actual_count = count(age -> age_group_idx(age, event_config.age_boundaries) == i, ages)
                actual = n > 0 ? actual_count / n : 0.0
                println("      <=$boundary: expected $(fmt(expected*100))% actual $(fmt(actual*100))% ($actual_count/$n)")
            end
            i = length(event_config.age_boundaries) + 1
            expected = expected_age[i]
            actual_count = count(age -> age_group_idx(age, event_config.age_boundaries) == i, ages)
            actual = n > 0 ? actual_count / n : 0.0
            println("      >$(event_config.age_boundaries[end]): expected $(fmt(expected*100))% actual $(fmt(actual*100))% ($actual_count/$n)")

            println("    Sex distribution:")
            for (s, label) in [(1, "Male"), (2, "Female")]
                expected = expected_sex[s]
                actual_count = count(==(s), sexes)
                actual = n > 0 ? actual_count / n : 0.0
                println("      $label: expected $(fmt(expected*100))% actual $(fmt(actual*100))% ($actual_count/$n)")
            end
        end
    end
end


# Compares each event's actual compartment counts against an analytically
# expected value (population-wide prevalence scaled to event size, plus a
# transmission-probability estimate for infected_at_event), as a Z-score.
# The transmission rate used is the mean transmission_prob across this
# event's actual attendees, not a flat population-wide constant — this
# matters when categories concentrate superspreaders (via min_superspreaders
# or core targeting), since a flat rate would systematically underestimate
# expected infections for those events.
#
# The compartment baseline assumes attendees are a random sample of the
# population. Core and loyal repeat attendees had extra exposure at earlier
# draws, so in later draws recovered/infectious can sit above the baseline:
# large Z-scores there reflect real interaction, not an error.
#
# The Infected row's std is binomial, sqrt(S*p*(1-p)), so it has a Z-score
# even with one replicate. It is a lower bound: it ignores variation in the
# number of infectious attendees, contact overdispersion and superspreaders,
# which inflates Z. With several replicates, observed is a mean but the std
# is still per-run, which deflates Z.
function validate_epidemic_state(aggregated, events, bd, event_config, people, attendees)
    println("\n=== Epidemic State Validation ===")

    pop_size = nrow(people)
    unique_dates = unique(e.date for e in events)

    counts_per_date = Dict(date => NamedTuple[] for date in unique_dates)
    for rd in runs(bd)
        inf_log = infections(rd)
        for date in unique_dates
            push!(counts_per_date[date], population_compartment_counts(inf_log, date, pop_size))
        end
    end

    pop_counts = Dict(
        date => (
            infectious=mean(c.infectious for c in cs),
            exposed=mean(c.exposed for c in cs),
            recovered=mean(c.recovered for c in cs),
            dead=mean(c.dead for c in cs),
            same_day_other=mean(c.same_day_other for c in cs),
        )
        for (date, cs) in counts_per_date
    )

    # z_str is "n/a" rather than a Z-score of 0 when std == 0 (e.g. with a
    # single simulation replicate), since 0 would misleadingly look like a
    # perfect match rather than "no variance to compare against."
    function print_validation_row(label, expected, observed, std)
        z_str = std > 0 ? fmt((observed - expected) / std) : "n/a"
        println("  $(rpad(label, 20)) $(rpad(fmt(expected), 12)) $(rpad(fmt(observed), 12)) $z_str")
    end

    for event in events
        metrics = aggregated[event.id]
        pc = pop_counts[event.date]

        expected_infectious = (pc.infectious / pop_size) * event.n
        expected_exposed = (pc.exposed / pop_size) * event.n
        expected_recovered = (pc.recovered / pop_size) * event.n
        expected_dead = (pc.dead / pop_size) * event.n
        # observed susceptible also excludes same-day infections in other
        # settings (GlobalSetting runs last), so the baseline must too
        expected_susceptible = event.n - expected_infectious - expected_exposed -
                               expected_recovered - expected_dead -
                               (pc.same_day_other / pop_size) * event.n

        n_section = event.n
        event_transmission_rate = mean(people.transmission_prob[idx] for idx in attendees[event.id])
        exponent = n_section > 1 ?
                   metrics[:infectious].mean * event.mean_contacts *
                   event_transmission_rate / (n_section - 1) : 0.0
        # cross-section contacts: an infectious person in section t draws them
        # from everyone in the same draw outside t, i.e. (draw total - n_t) people
        if event.cross_section_mean_contacts > 0
            siblings = filter(e -> e.category_id == event.category_id && e.draw_id == event.draw_id &&
                                   e.section_id != event.section_id, events)
            n_draw = event.n + sum(e.n for e in siblings; init=0)
            for e in siblings
                sibling_rate = mean(people.transmission_prob[idx] for idx in attendees[e.id])
                exponent += aggregated[e.id][:infectious].mean * event.cross_section_mean_contacts *
                            sibling_rate / (n_draw - e.n)
            end
        end
        p_infected = 1 - exp(-exponent)
        expected_infected = metrics[:susceptible].mean * p_infected
        std_infected = sqrt(metrics[:susceptible].mean * p_infected * (1 - p_infected))

        println("\nEvent $(event.id) — Day $(event.date) — n=$(event.n)")
        println("  $(rpad("", 20)) $(rpad("Expected", 12)) $(rpad("Observed", 12)) Z-score")
        println("  " * "-"^55)
        print_validation_row("Susceptible", expected_susceptible, metrics[:susceptible].mean, metrics[:susceptible].std)
        print_validation_row("Infectious", expected_infectious, metrics[:infectious].mean, metrics[:infectious].std)
        print_validation_row("Exposed", expected_exposed, metrics[:exposed].mean, metrics[:exposed].std)
        print_validation_row("Recovered", expected_recovered, metrics[:recovered].mean, metrics[:recovered].std)
        print_validation_row("Dead", expected_dead, metrics[:dead].mean, metrics[:dead].std)
        println("  " * "-"^55)
        print_validation_row("Infected", expected_infected, metrics[:infected_at_event].mean, std_infected)
    end
end


# === Auto-run validation ===

validate_assignment(people, events, event_config, attendees)
validate_category_attendance(events, event_config, attendees, core_groups)
validate_results(aggregated, events)
validate_demographics(people, events, event_config, attendees)
validate_epidemic_state(aggregated, events, bd, event_config, people, attendees)

println("End Validation")