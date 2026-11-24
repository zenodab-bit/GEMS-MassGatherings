# ===========================================================================
# contacts.jl
#
# Extends GEMS's sample_contacts! for two ContactSamplingMethods, dispatched
# by the "type" set in each setting's TOML config: NegBinContacts for
# ordinary settings (household, office, etc.), EventContacts specifically
# for GlobalSetting (mass gatherings).
# ===========================================================================

import GEMS: sample_contacts!


# === NegBin Contact Sampling for normal settings ===

@with_kw mutable struct NegBinContacts <: GEMS.ContactSamplingMethod
    mean_contacts::Float64
    std_contacts::Float64
end

# Draws n contacts (via sample_n_contacts) from present_individuals,
# excluding ego, using distinct offsets so no contact repeats and none is ego.
function GEMS.sample_contacts!(
    indivs::Vector{Individual},
    csm::NegBinContacts,
    setting::Setting,
    individual_index::Int,
    present_individuals::Vector{Individual},
    tick::Int16,
    replace::Bool,
    rng::Xoshiro
)
    empty!(indivs)
    length(present_individuals) <= 1 && return indivs

    n_contacts = min(sample_n_contacts(rng, csm.mean_contacts, csm.std_contacts), length(present_individuals) - 1)
    n_contacts <= 0 && return indivs

    # distinct offsets => distinct contacts, and offset never lands on ego (offset ∈ 1:length-1, mod length)
    offsets = sample(rng, 1:(length(present_individuals)-1), n_contacts, replace=false)

    resize!(indivs, n_contacts)
    for i in 1:n_contacts
        contact_index = mod(individual_index + offsets[i] - 1, length(present_individuals)) + 1
        indivs[i] = present_individuals[contact_index]
    end

    return indivs
end


# === Event Contact Sampling ===

# rosters is rebuilt once per tick and pass (see below), shared across every
# call of that pass rather than rescanning the whole population per call.
#
# pass: which GlobalSetting pass of the day is running. Pass 1 is GEMS's own
# (normal events, and the first stage of staged ones); pass s is stage s of
# staged categories, run by daily_global_setting! at the end of the day.
# pass_infections records, for passes after the first, who was infected in
# which pass: (person id, day) => pass. Infections at events not in it
# happened in pass 1.
@with_kw mutable struct EventContacts <: GEMS.ContactSamplingMethod
    cached_tick::Int16 = Int16(-1)
    cached_pass::Int = 0
    pass::Int = 1
    rosters::Dict{NTuple{4,Int32},Vector{Individual}} = Dict{NTuple{4,Int32},Vector{Individual}}()
    draw_rosters::Dict{NTuple{3,Int32},Vector{Individual}} = Dict{NTuple{3,Int32},Vector{Individual}}()
    pass_infections::Dict{Tuple{Int32,Int16},Int} = Dict{Tuple{Int32,Int16},Int}()
end

# The GlobalSetting pass an attendance belongs to: stage s of a staged
# category in pass s, normal events (stage 0) in pass 1.
event_pass(stage_id) = max(Int(stage_id), 1)

# The position of x's attendance on day `tick` in pass `pass`, or nothing.
event_entry(x::Individual, tick::Int16, pass::Int) =
    findfirst(k -> x.event_dates[k] == tick && event_pass(x.stage_ids[k]) == pass, eachindex(x.event_dates))

# Attendance is fixed before the simulation, but on the day itself people
# who are dead, severely ill or in hospital don't go: they're left out of
# that day's rosters, so they neither infect nor take up contact slots.
# (GEMS already stops them infecting outside the household, but would still
# let others draw them as contacts.) People with milder symptoms still go.
# Uses GEMS's state as of `tick`, computed from each person's disease dates.
absent_from_events(x::Individual, tick::Int16) =
    GEMS.is_dead(x, tick) || GEMS.is_severe(x, tick) || GEMS.is_hospitalized(x, tick)

# Finds ego's event (if any) for this tick and pass, then draws contacts from
# others in the same (category, draw, stage, section); for a staged
# category the section is the venue. Builds/caches a roster of every
# attendee keyed that way the first time this tick and pass are seen, so
# repeated calls reuse it instead of rescanning present_individuals each time.
function GEMS.sample_contacts!(
    indivs::Vector{Individual},
    event_contacts::EventContacts,
    setting::GEMS.GlobalSetting,
    individual_index::Int,
    present_individuals::Vector{Individual},
    tick::Int16,
    replace::Bool,
    rng::Xoshiro
)
    empty!(indivs)
    ego = present_individuals[individual_index]

    pass = event_contacts.pass
    if event_contacts.cached_tick != tick || event_contacts.cached_pass != pass
        empty!(event_contacts.rosters)
        empty!(event_contacts.draw_rosters)
        n_absent = 0
        for x in present_individuals
            isempty(x.event_dates) && continue
            jdx = event_entry(x, tick, pass)
            jdx === nothing && continue
            if absent_from_events(x, tick)
                n_absent += 1
                continue
            end
            section_key = (x.category_ids[jdx], x.draw_ids[jdx], x.stage_ids[jdx], x.section_ids[jdx])
            draw_key = (x.category_ids[jdx], x.draw_ids[jdx], x.stage_ids[jdx])
            push!(get!(() -> Individual[], event_contacts.rosters, section_key), x)
            push!(get!(() -> Individual[], event_contacts.draw_rosters, draw_key), x)
        end
        event_contacts.cached_tick = tick
        event_contacts.cached_pass = pass
        n_absent > 0 && @debug "Day $tick, pass $pass: $n_absent assigned attendee(s) absent (dead, severe or hospitalised)"
    end

    idx = event_entry(ego, tick, pass)
    idx === nothing && return indivs
    absent_from_events(ego, tick) && return indivs

    section_key = (ego.category_ids[idx], ego.draw_ids[idx], ego.stage_ids[idx], ego.section_ids[idx])
    draw_key = (ego.category_ids[idx], ego.draw_ids[idx], ego.stage_ids[idx])

    section_roster = get(event_contacts.rosters, section_key, nothing)
    same_section_others = section_roster === nothing ? Individual[] : filter(x -> x !== ego, section_roster)

    # within-section contacts, using this attendee's own section rate
    n_within = min(sample_n_contacts(rng, ego.mean_event_contacts[idx], ego.std_event_contacts[idx]), length(same_section_others))
    if n_within > 0
        append!(indivs, sample(rng, same_section_others, n_within, replace=false))
    end

    # cross-section contacts, using the category's cross-section rate (0 by
    # default keeps sections fully isolated, matching prior behavior)
    if ego.cross_section_mean_contacts[idx] > 0
        draw_roster = get(event_contacts.draw_rosters, draw_key, nothing)
        if draw_roster !== nothing
            same_section_set = Set(section_roster === nothing ? Individual[] : section_roster)
            cross_others = filter(x -> x !== ego && x ∉ same_section_set, draw_roster)
            n_cross = min(sample_n_contacts(rng, ego.cross_section_mean_contacts[idx], ego.cross_section_std_contacts[idx]), length(cross_others))
            if n_cross > 0
                append!(indivs, sample(rng, cross_others, n_cross, replace=false))
            end
        end
    end

    return indivs
end


# === GlobalSetting membership per day ===

# GEMS puts the whole population in the GlobalSetting, but on any given day
# only that day's event attendees do anything there (everyone else gets no
# event contacts and transmission 0). So the GlobalSetting is given only
# that day's attendees, and nobody on days without events: same results,
# much less work per day.
#
# Passed to Batch as `setup`, so it runs once per simulation, before day 0.
# It builds one list per event day (all events and sections of that day
# together; EventContacts splits them by category, draw and section), then
# installs a stepmod that GEMS runs at the end of every day to put the next
# day's list in place. Lists keep population order, so the same seed gives
# the same results as with the whole population in the GlobalSetting.
#
# expected_attendances is the total number of event attendances in the
# population table (sum of length(event_dates)). If GEMS's individuals hold
# a different total, the event columns didn't reach GEMS intact; that would
# otherwise fail silently (no event contacts at all), so it stops with an error.
#
# events_on = false leaves the GlobalSetting empty on every day, so no event
# takes place, while everything else (events, attendees, transmission_prob)
# stays as set up. That's the run without events that run_comparison
# (comparison.jl) compares against. Before the first event day the
# GlobalSetting is empty either way, so with the same seeds both runs are
# identical up to that day.
#
# Staged categories: each stage gets its own pass through the GlobalSetting.
# GEMS's own pass of the day holds pass 1 (normal events and the first
# stage); at the end of the day (GEMS runs the stepmod before moving to the
# next day, so it's still the same day), each further stage s is put in the
# GlobalSetting and GEMS's spread_infection! runs it once more. Who gets
# infected in those passes is recorded in the EventContacts' pass_infections.
# The order of passes doesn't matter for who can infect whom (people infected
# today only become infectious days later); someone infected in an earlier
# pass just can't be infected again. GEMS never switches the GlobalSetting
# off (deactivate! is a no-op for it), so the extra passes don't affect the
# next day. Without staged categories there's only pass 1, as before.
function daily_global_setting!(sim, expected_attendances::Integer; events_on::Bool = true)
    gs = settings(sim, GlobalSetting)[1]
    csm = gs.contact_sampling_method

    # (day, pass) => that day's attendees in that pass; people without events
    # never enter a list
    by_day = Dict{Tuple{Int16,Int},Vector{Individual}}()
    for ind in gs.individuals
        for (d, s) in zip(ind.event_dates, ind.stage_ids)
            push!(get!(() -> Individual[], by_day, (Int16(d), event_pass(s))), ind)
        end
    end
    nobody = Individual[]
    last_pass = Dict{Int16,Int}()   # day => its number of passes
    for (d, p) in keys(by_day)
        last_pass[d] = max(get(last_pass, d, 1), p)
    end

    found = sum(length, values(by_day); init=0)
    found == expected_attendances ||
        error("GEMS's individuals hold $found event attendances, expected $expected_attendances: " *
              "the event columns did not reach GEMS intact.")

    set_list!(t, p) = (gs.individuals = events_on ? get(by_day, (t, p), nobody) : nobody)

    # the passes after GEMS's own, on day t
    function extra_passes!(s, t)
        (events_on && !GEMS.is_dormant(s)) || return
        for p in 2:get(last_pass, t, 1)
            set_list!(t, p)
            csm.pass = p
            not_yet = [ind for ind in gs.individuals if !GEMS.infected(ind)]
            GEMS.spread_infection!(gs, s, GEMS.pathogen(s))
            for ind in not_yet
                GEMS.infected(ind) && (csm.pass_infections[(Int32(GEMS.id(ind)), t)] = p)
            end
        end
        csm.pass = 1
    end

    set_list!(tick(sim), 1)                # first day
    previous = sim.stepmod                 # keep any stepmod already set
    sim.stepmod = s -> (previous(s); extra_passes!(s, tick(s)); set_list!(tick(s) + Int16(1), 1))
    return sim
end

