# ===========================================================================
# validation.jl
#
## Checks on a finished run. Including this file only defines functions;
# run_validation_checks(...) at the bottom runs them all.
# ===========================================================================

# Population-wide SEIR compartment counts as of `tick`, computed directly
# from the infection log (not from tick_cases, which turned out to report
# same-day incidence rather than cumulative prevalence). Used by
# validate_epidemic_state as the population-wide baseline. Uses the same
# compartments as analyze_event_population: absent (dead, severe or
# hospitalised that day) first, then infectious / recovered / exposed.
function population_compartment_counts(inf_log::DataFrame, tick::Int, pop_size::Int)
    infectious = 0
    exposed = 0
    recovered = 0
    absent = 0
    same_day_other = 0   # infected on `tick` outside GlobalSetting, before events run

    for row in eachrow(inf_log)
        if row.tick == tick
            row.setting_type != 'g' && (same_day_other += 1)
            continue
        end
        row.tick > tick && continue
        if absent_on(tick, row)
            absent += 1
        elseif row.infectiousness_onset <= tick &&
           (row.recovery > tick || row.recovery == -1)
            infectious += 1
        elseif row.recovery != -1 && row.recovery <= tick
            recovered += 1
        else
            exposed += 1
        end
    end

    susceptible = pop_size - infectious - exposed - recovered - absent - same_day_other
    return (susceptible=susceptible, infectious=infectious, exposed=exposed,
        recovered=recovered, absent=absent, same_day_other=same_day_other)
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

    # nobody may attend two events on the same day, except different stages
    # of one draw of a staged category
    n_double = count(1:nrow(people)) do i
        ds = people.event_dates[i]
        allunique(ds) && return false
        any(unique(ds)) do d
            ks = findall(==(d), ds)
            length(ks) > 1 && !(all(k -> people.stage_ids[i][k] > 0, ks) && allunique(people.stage_ids[i][ks]) &&
                                allequal(people.category_ids[i][ks]) && allequal(people.draw_ids[i][ks]))
        end
    end
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
        if category isa StagedCategory
            validate_stage_attendance(category, cat_events, attendees)
            continue
        end

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

# A staged category's attendance, per draw and stage: attendees (all venues
# together) against the stage's size, and the overlap reached against the
# one set: how many of the stage's attendees were at the previous stage of
# the same draw, wanted round(overlap × size).
function validate_stage_attendance(category::StagedCategory, cat_events, attendees)
    println("\n  $(category.name) (id=$(category.id), staged, n_draws=$(category.n_draws)):")
    for draw in sort(unique(e.draw_id for e in cat_events))
        previous = Set{Int}()
        for (s, stage) in enumerate(category.stages)
            evs = filter(e -> e.draw_id == draw && e.stage_id == s, cat_events)
            isempty(evs) && continue
            members = reduce(vcat, [attendees[e.id] for e in evs])
            size_ = sum(e.n for e in evs)
            label = "    draw $draw (day $(evs[1].date)), stage $s ($(stage.name)): $(length(members)) of $size_ places"
            if s == 1
                println(label)
            else
                wanted = round(Int, stage.overlap * size_)
                got = count(in(previous), members)
                println(label, "; overlap $got from stage $(s - 1), wanted $wanted ", got == wanted ? "✓" : "✗")
            end
            allunique(members) || println("    ✗ someone is at two venues of this stage")
            previous = Set(members)
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

        # susceptible + infectious + exposed + recovered + absent + same-day-other should equal section size
        total = metrics[:susceptible].mean + metrics[:infectious].mean +
                metrics[:exposed].mean + metrics[:recovered].mean +
                metrics[:absent].mean + metrics[:same_day_other].mean
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

        # a staged category's stages each have their own weights: one group per stage
        groups = category isa StagedCategory ?
            [(stage, "\nStaged category $(category.id), stage $s ($(stage.name)):",
              filter(e -> e.stage_id == s, cat_events)) for (s, stage) in enumerate(category.stages)] :
            [(category, "\nCategory $(category.id) — loyalty=$(category.loyalty):", cat_events)]

        for (weights_from, header, group_events) in groups
            # expected share per (age group, sex) cell, built the same way
            # compute_weights builds sampling weights: age weight x sex weight
            n_groups = length(event_config.age_boundaries) + 1
            mass = [(isempty(weights_from.age_weights) ? event_config.age_dist[g] : weights_from.age_weights[g]) *
                    (isempty(weights_from.sex_weights) ? event_config.sex_dist[g][s] : weights_from.sex_weights[s])
                    for g in 1:n_groups, s in 1:2]
            mass ./= sum(mass)
            expected_age = vec(sum(mass, dims=2))
            expected_sex = vec(sum(mass, dims=1))

            println(header)

            for event in group_events
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
end


# --- Expected infections at events (the Infected row) ---
#
# A check of the event mechanism (contacts.jl, transmission.jl) that
# doesn't reuse its code: the simulation draws contacts and infections at
# random; this computes, for each replicate, the probability of each
# outcome from the rules the events are meant to follow, given who was
# present and infectious on the event day:
#   - an infectious attendee draws K contacts (negative binomial with the
#     section's mean and std, Poisson when std^2 <= mean), capped at the
#     number of others present in their section, and meets that many
#     different people there, chosen at random;
#   - with cross-section contacts, the same again with the category's
#     cross-section mean and std, among those present in the draw's other
#     sections;
#   - each contact infects a susceptible with the infecter's own
#     transmission_prob.
# If the simulation doesn't follow these rules (e.g. contacts from the
# wrong section, a wrong number of contacts, the wrong person's
# transmission_prob), the observed infections drift away from these.
# Inputs both share (the attendee lists, transmission_prob, the configured
# contact numbers) aren't tested here.

# Expected contacts of one person, E[min(K, m)]: K drawn with the given
# mean and std, capped at m, the number of people there are to meet.
# Written from the distributions' definitions, independently of
# sample_n_contacts / negbin_params (utils.jl), so a mistake there
# can't hide in the check. A negative binomial with mean μ and variance
# σ² > μ has success probability μ/σ² and r = μ²/(σ² - μ).
function expected_capped_contacts(mean_k, std_k, m)
    (mean_k <= 0 || m <= 0) && return 0.0
    v = std_k^2
    K = v <= mean_k ? Poisson(mean_k) : NegativeBinomial(mean_k^2 / (v - mean_k), mean_k / v)
    return sum(ccdf(K, k) for k in 0:(m - 1))   # E[min(K, m)] = Σ P(K > k), k < m
end

# The state on day d of the people in `ids`, from one replicate's infection
# log, using each person's latest infection up to day d: :absent (dead,
# severe or hospitalised that day, see absent_on), :infectious (infected
# before d, infectious on d and not yet recovered), or :other (infected but
# neither). People not in the result were never infected up to day d.
function attendee_states(inf_log::DataFrame, d::Integer, ids::Set{Int32})
    latest = Dict{Int32,Int}()   # id => row of the latest infection up to day d
    for (i, (t, b)) in enumerate(zip(inf_log.tick, inf_log.id_b))
        (t <= d && b in ids) || continue
        j = get(latest, b, 0)
        (j == 0 || inf_log.tick[j] <= t) && (latest[b] = i)
    end
    states = Dict{Int32,Symbol}()
    for (b, i) in latest
        row = inf_log[i, :]
        states[b] = absent_on(d, row) ? :absent :
                    (row.tick < d && 0 <= row.infectiousness_onset <= d &&
                     (row.recovery > d || row.recovery == -1)) ? :infectious : :other
    end
    return states
end

# Per event and replicate: the expected number infected at the event. Per
# replicate r, a susceptible attendee present in section s escapes infection
# with probability
#   Π over infectious i in s:            1 - c_i / m_s * p_i
#   Π over infectious i in another section t of the draw (cross-section):
#                                        1 - c'_i / m'_t * p_i
# (c = expected capped contacts, m = how many others they can meet, p_i = i's
# transmission_prob); expected infected = S_r * (1 - escape), S_r the
# susceptible attendees present (from event_results). Returns event id =>
# one expected value per replicate.
function expected_infected_at_events(events, people, setups, bd, event_results)
    expected = Dict(e.id => Float64[] for e in events)
    # one group per draw, and for staged categories per stage of it: each
    # stage is its own GlobalSetting pass, with contacts only within the venue
    draws = unique((e.category_id, e.draw_id, e.stage_id) for e in events)

    for (r, rd) in enumerate(runs(bd))
        s = setups[length(setups) == 1 ? 1 : r]
        inf_log = infections(rd)
        for (cat, draw, stage) in draws
            group = filter(e -> e.category_id == cat && e.draw_id == draw && e.stage_id == stage, events)
            d = group[1].date
            ids = Set{Int32}(people.id[idx] for e in group for idx in s.attendees[e.id])
            states = attendee_states(inf_log, d, ids)
            state(idx) = get(states, Int32(people.id[idx]), :never)
            present = Dict(e.id => [idx for idx in s.attendees[e.id] if state(idx) != :absent] for e in group)
            infectious = Dict(e.id => [idx for idx in present[e.id] if state(idx) == :infectious] for e in group)
            n_draw = sum(length(present[e.id]) for e in group)

            for e in group   # the susceptible's section
                escape = 1.0
                m = length(present[e.id]) - 1
                if m > 0
                    c = expected_capped_contacts(e.mean_contacts, e.std_contacts, m)
                    for i in infectious[e.id]
                        escape *= 1 - c / m * s.transmission_prob[i]
                    end
                end
                for t in group   # infecters in the draw's other sections
                    (t.id != e.id && t.cross_section_mean_contacts > 0) || continue
                    m2 = n_draw - length(present[t.id])
                    m2 > 0 || continue
                    c2 = expected_capped_contacts(t.cross_section_mean_contacts, t.cross_section_std_contacts, m2)
                    for i in infectious[t.id]
                        escape *= 1 - c2 / m2 * s.transmission_prob[i]
                    end
                end
                push!(expected[e.id], event_results[r][e.id].susceptible * (1 - escape))
            end
        end
    end
    return expected
end

# Compares each event's actual compartment counts against an expected
# value, as a Z-score.
#
# Susceptible to Absent: population-wide prevalence scaled to event size.
# This assumes attendees are a random sample of the population. Core and
# loyal repeat attendees had extra exposure at earlier draws, so in later
# draws recovered/infectious can sit above the baseline: large Z-scores
# there reflect real interaction, not an error. Observed is the mean over
# replicates, so each of these rows' std (one run's spread) is divided by
# sqrt(n_runs) to get the spread of that mean. Read them as sanity checks.
#
# Infected: the expected infections at the event, computed per replicate
# from who was present and infectious (expected_infected_at_events above).
# The Z-score uses the per-replicate differences, observed minus expected:
# their mean over their standard error (std / √replicates). The spread
# comes from the replicates themselves, since infections at an event are
# far more variable than binomial (one infectious person with many contacts
# can infect several at once); so it needs at least 2 replicates, and more
# give a sharper test. This one is a real test of the event code: |Z| above
# about 3 is worth a closer look.
function validate_epidemic_state(aggregated, events, bd, event_config, people, setups, event_results)
    println("\n=== Epidemic State Validation ===")

    n_runs = length(runs(bd))
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
            absent=mean(c.absent for c in cs),
            same_day_other=mean(c.same_day_other for c in cs),
        )
        for (date, cs) in counts_per_date
    )

    # z_str is "n/a" rather than a Z-score of 0 when std == 0 (e.g. with a
    # single simulation replicate), since 0 would misleadingly look like a
    # perfect match rather than "no variance to compare against."
    # std is one run's spread; std / sqrt(n_runs) is the spread of the mean.
    function print_validation_row(label, expected, observed, std)
        se = std / sqrt(n_runs)
        z_str = se > 0 ? fmt((observed - expected) / se) : "n/a"
        println("  $(rpad(label, 20)) $(rpad(fmt(expected), 12)) $(rpad(fmt(observed), 12)) $z_str")
    end

    infected_expected = expected_infected_at_events(events, people, setups, bd, event_results)

    for event in events
        metrics = aggregated[event.id]
        pc = pop_counts[event.date]

        expected_infectious = (pc.infectious / pop_size) * event.n
        expected_exposed = (pc.exposed / pop_size) * event.n
        expected_recovered = (pc.recovered / pop_size) * event.n
        expected_absent = (pc.absent / pop_size) * event.n
        # observed susceptible also excludes same-day infections in other
        # settings (GlobalSetting runs last), so the baseline must too
        expected_susceptible = event.n - expected_infectious - expected_exposed -
                               expected_recovered - expected_absent -
                               (pc.same_day_other / pop_size) * event.n

        println("\nEvent $(event.id) — Day $(event.date) — n=$(event.n)")
        println("  $(rpad("", 20)) $(rpad("Expected", 12)) $(rpad("Observed", 12)) Z-score")
        println("  " * "-"^55)
        print_validation_row("Susceptible", expected_susceptible, metrics[:susceptible].mean, metrics[:susceptible].std)
        print_validation_row("Infectious", expected_infectious, metrics[:infectious].mean, metrics[:infectious].std)
        print_validation_row("Exposed", expected_exposed, metrics[:exposed].mean, metrics[:exposed].std)
        print_validation_row("Recovered", expected_recovered, metrics[:recovered].mean, metrics[:recovered].std)
        print_validation_row("Absent", expected_absent, metrics[:absent].mean, metrics[:absent].std)
        println("  " * "-"^55)
        # per replicate: observed minus expected; Z = their mean / standard error
        expected_r = infected_expected[event.id]
        diffs = [event_results[r][event.id].infected_at_event - expected_r[r] for r in eachindex(expected_r)]
        se = length(diffs) > 1 ? std(diffs) / sqrt(length(diffs)) : 0.0
        z_str = se > 0 ? fmt(mean(diffs) / se) : "n/a"
        println("  $(rpad("Infected", 20)) $(rpad(fmt(mean(expected_r)), 12)) " *
                "$(rpad(fmt(metrics[:infected_at_event].mean), 12)) $z_str")
    end
end


# === Run everything ===

# Runs every check, writes the output to validation.txt in run_folder, then
# shows it on screen too.
function run_validation_checks(; people, events, event_config, attendees,
                                 core_groups, aggregated, event_results, bd, run_folder,
                                 setups = [(; attendees, core_groups, transmission_prob = people.transmission_prob)])
    validation_file = joinpath(run_folder, "validation.txt")
    open(validation_file, "w") do io
        redirect_stdout(io) do
            length(setups) > 1 &&
                println("Each replicate drew its own setup (vary_per_replicate). The assignment, attendance " *
                        "and demographics checks below use replicate 1's; every replicate's is drawn with " *
                        "the same rules. The epidemic-state check uses all of them.")
            validate_assignment(people, events, event_config, attendees)
            validate_category_attendance(events, event_config, attendees, core_groups)
            validate_results(aggregated, events)
            validate_demographics(people, events, event_config, attendees)
            validate_epidemic_state(aggregated, events, bd, event_config, people, setups, event_results)
        end
    end
    print(read(validation_file, String))
    println("Saved: validation")
end

