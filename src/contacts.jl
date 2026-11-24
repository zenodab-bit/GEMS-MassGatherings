# ===========================================================================
# contacts.jl
#
# GEMS's sample_contacts! for the two contact methods named in the TOML:
# NegBinContacts for normal settings, EventContacts for the GlobalSetting
# (events). Also daily_global_setting!, which puts each day's attendees in
# the GlobalSetting and runs the extra passes for staged categories.
# ===========================================================================

import GEMS: sample_contacts!


# === NegBin Contact Sampling for normal settings ===

@with_kw mutable struct NegBinContacts <: GEMS.ContactSamplingMethod
    mean_contacts::Float64
    std_contacts::Float64
end

# n distinct contacts (n from sample_n_contacts) among the others present:
# distinct offsets from ego, so no repeats and never ego itself.
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

    offsets = sample(rng, 1:(length(present_individuals)-1), n_contacts, replace=false)

    resize!(indivs, n_contacts)
    for i in 1:n_contacts
        contact_index = mod(individual_index + offsets[i] - 1, length(present_individuals)) + 1
        indivs[i] = present_individuals[contact_index]
    end

    return indivs
end


# === Event Contact Sampling ===

# rosters: attendees present by (category, draw, stage, section), built
# once per day and pass. pass: the running pass (stage s in pass s, normal
# events in 1). pass_infections: (person id, day) => pass, for passes > 1.
@with_kw mutable struct EventContacts <: GEMS.ContactSamplingMethod
    cached_tick::Int16 = Int16(-1)
    cached_pass::Int = 0
    pass::Int = 1
    rosters::Dict{NTuple{4,Int32},Vector{Individual}} = Dict{NTuple{4,Int32},Vector{Individual}}()
    draw_rosters::Dict{NTuple{3,Int32},Vector{Individual}} = Dict{NTuple{3,Int32},Vector{Individual}}()
    pass_infections::Dict{Tuple{Int32,Int16},Int} = Dict{Tuple{Int32,Int16},Int}()
end

# The pass an attendance belongs to: stage s in pass s, normal events in 1.
event_pass(stage_id) = max(Int(stage_id), 1)

# Position of x's attendance on day `tick` in pass `pass`, or nothing.
event_entry(x::Individual, tick::Int16, pass::Int) =
    findfirst(k -> x.event_dates[k] == tick && event_pass(x.stage_ids[k]) == pass, eachindex(x.event_dates))

# Dead, severely ill or hospitalised people don't go to their events that
# day: they neither infect nor take up contact slots (GEMS alone would
# still let others draw them as contacts). Milder cases still go.
absent_from_events(x::Individual, tick::Int16) =
    GEMS.is_dead(x, tick) || GEMS.is_severe(x, tick) || GEMS.is_hospitalized(x, tick)

# Contacts of ego at its event this day and pass: within its section (or
# venue) at its own rate, plus cross-section contacts within the draw.
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

    n_within = min(sample_n_contacts(rng, ego.mean_event_contacts[idx], ego.std_event_contacts[idx]), length(same_section_others))
    if n_within > 0
        append!(indivs, sample(rng, same_section_others, n_within, replace=false))
    end

    # cross-section contacts at the category's rate (0 keeps sections apart)
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

# The GlobalSetting gets only each day's attendees: same results, much less
# work. Run once per simulation as the Batch's `setup`: builds a list per
# (day, pass) and a stepmod that swaps in the next day's list at the end of
# each day. The lists keep population order, which is what makes the
# results identical to having everyone in the GlobalSetting. Errors if the
# event columns didn't reach GEMS intact. events_on = false keeps it empty
# (the run without events).
#
# Stage s > 1 of a staged category runs as an extra pass in the stepmod,
# with its infections recorded in pass_infections. Pass order doesn't
# matter: people infected today only become infectious later. GEMS's
# deactivate! does nothing for the GlobalSetting, so the extra passes
# can't switch it off for the next day.
function daily_global_setting!(sim, expected_attendances::Integer; events_on::Bool = true)
    gs = settings(sim, GlobalSetting)[1]
    csm = gs.contact_sampling_method

    by_day = Dict{Tuple{Int16,Int},Vector{Individual}}()   # (day, pass) => attendees
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

