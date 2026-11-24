# ===========================================================================
# MassGatherings.jl
#
# The package: mass gathering events on top of the GEMS epidemic simulator.
# Files load in this order (each may use the earlier ones). The exports are
# what users call; the rest is reachable as MassGatherings.name.
# ===========================================================================

module MassGatherings

using GEMS, Parameters, DataFrames, Distributions, CSV, CategoricalArrays,
    JLD2, Random, StatsBase, Statistics, Plots, StatsPlots, Dates, TOML

include("utils.jl")          # distributions, age groups, plot colours, tables
include("config.jl")         # configuration structs and their checks
include("events.jl")         # drawing the events
include("attendance.jl")     # preparing the population, assigning attendees
include("contacts.jl")       # contacts, daily GlobalSetting lists and passes
include("transmission.jl")   # transmission
include("results.jl")        # per-event counts, aggregation, saving
include("scenario.jl")       # run_scenario
include("comparison.jl")     # with and without events
include("sweeps.jl")         # sweeps over one setting
include("analysis.jl")       # text metrics and plots
include("validation.jl")     # checks of a finished run

# configuring events (config/user_config.jl)
export EventSection, Category, Venue, Stage, StagedCategory, EventConfig, validate_config
export gamma_params, beta_params
# running (run.jl)
export run_scenario, run_comparison, run_sweep, with_date, with_date_range
# GEMS finds these by the names written in the TOML, which only match
# without the "MassGatherings." prefix, i.e. when exported
export NegBinContacts, EventContacts, SettingRate

end # module
