# ===========================================================================
# 5_Transmission.jl
#
# Extends GEMS's transmission_probability for SettingRate: in every setting,
# events included, a contact infects with the infecter's own
# transmission_prob, unless the contact has already recovered (natural
# immunity: without this check GEMS would let recovered people be infected
# again).
#
# Who meets whom at events is decided entirely in 4_Contacts.jl: contacts
# come only from the attendee's own section (plus cross-section contacts
# within the same draw), so there's no "same event" check here. The exact
# expected-infections check in 7_Validation.jl tests that this holds.
# ===========================================================================

import GEMS.transmission_probability

# Transmission itself uses each person's transmission_prob (set in
# 3_Population.jl), so there's no rate parameter here. GEMS first builds one
# from the TOML (its `parameters` table must exist but can be empty), then
# run_scenario (0_Helpers.jl) passes its own.
struct SettingRate <: GEMS.TransmissionFunction end

function GEMS.transmission_probability(
    transFunc::SettingRate,
    infecter::Individual,
    infected::Individual,
    setting::Setting,
    tick::Int16
)::Float64

    if -1 < recovery(infected) <= tick
        return 0.0
    end

    return infecter.transmission_prob
end


println("End Transmission")
