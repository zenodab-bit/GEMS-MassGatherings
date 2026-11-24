# ===========================================================================
# transmission.jl
#
# GEMS's transmission_probability for SettingRate: in every setting, a
# contact infects with the infecter's own transmission_prob, unless the
# contact has recovered (GEMS alone would allow reinfection). Who meets
# whom at events is decided in contacts.jl; the validation's Infected row
# checks it.
# ===========================================================================

import GEMS.transmission_probability

# No parameters: GEMS builds one from the TOML (whose `parameters` table
# must exist but can be empty), then run_scenario passes its own.
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

