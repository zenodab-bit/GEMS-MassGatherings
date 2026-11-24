# ===========================================================================
# MassGatherings.jl
#
# The package's module: mass gathering events (festivals, concerts, sport,
# staged events such as restaurants and a match) on top of the GEMS
# epidemic simulator. run.jl loads it with `using MassGatherings`.
#
# The files below are loaded in this order (each may use what the earlier
# ones define); the export list says what users are meant to call. Anything
# else is internal, reachable as MassGatherings.name if needed.
# ===========================================================================

module MassGatherings

using GEMS, Parameters, DataFrames, Distributions, CSV, CategoricalArrays,
    JLD2, Random, StatsBase, Statistics, Plots, StatsPlots, Dates, TOML

include("utils.jl")          # distributions, age groups, plot colours, table formatting
include("config.jl")         # the configuration structs and their checks
include("events.jl")         # drawing the events
include("attendance.jl")     # preparing the population, assigning attendees
include("contacts.jl")       # contacts at events and in other settings, daily passes
include("transmission.jl")   # transmission
include("results.jl")        # per-event counts, aggregation, saving
include("scenario.jl")       # run_scenario: one scenario from start to finish
include("comparison.jl")     # with and without events
include("sweeps.jl")         # sweeps over one setting
include("analysis.jl")       # text metrics and plots
include("validation.jl")     # checks of a finished run

# --- what users call ---
# configuring events (config/user_config.jl)
export EventSection, Category, Venue, Stage, StagedCategory, EventConfig, validate_config
export gamma_params, beta_params
# running (run.jl)
export run_scenario, run_comparison, run_sweep, with_date, with_date_range
# GEMS finds the contact and transmission types named in the TOML by their
# printed names; exported (and loaded with `using`), they print without the
# "MassGatherings." prefix and match "NegBinContacts", "EventContacts" and
# "SettingRate" as written there
export NegBinContacts, EventContacts, SettingRate

end # module
