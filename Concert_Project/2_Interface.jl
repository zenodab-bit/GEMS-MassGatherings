# ===========================================================================
# 2_Interface.jl
#
# Run this file to run a simulation. Loads packages, includes every other
# file, then calls run_scenario (0_Helpers.jl), which samples events and
# population from 1_UserConfig.jl's event_config, runs the GEMS batch, then
# triggers analysis and (optionally) validation.
#
# The settings below (n_simulations, run_validation, the comparison and
# sweep, and the population files) are the main things to adjust here;
# everything else you'd want to change lives in 1_UserConfig.jl.
# ===========================================================================

ENV["GKSwstype"] = "100"   # headless plotting backend; must be set before Plots loads

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GEMS, Parameters, DataFrames, Distributions, CSV,
    CategoricalArrays, JLD2, Random, StatsBase, Plots, 
    StatsPlots, Dates, TOML

# === Includes ===
include("0_Helpers.jl")
include("1_UserConfig.jl")
include("3_Population.jl")
include("4_Contacts.jl")
include("5_Transmission.jl")
include("6_Analysis.jl")
include("7_Validation.jl")

# === Run ===
n_simulations = 3
run_validation = true
# What each replicate draws for itself instead of sharing one draw (the
# events themselves always stay the same):
#   []                            - nothing: one setup for all replicates (default)
#   [:attendance]                 - who attends which event
#   [:transmission]               - everyone's transmission_prob
#   [:attendance, :transmission]  - both
vary_per_replicate = []
# true: run the scenario twice, with and without its events (same seeds),
# and compare the two (comparison.txt / .png). Takes twice as long.
compare_without_events = false

# Sweep: runs the comparison above once per value in sweep_values, to see
# how the events' effect depends on one setting (sweep.txt / .png). `sweep`
# says what changes for each value; anything it doesn't return stays as set
# here and in 1_UserConfig.jl. sweep_values = nothing: no sweep.
# Pick one of the two blocks below (or write your own):
sweep_values = nothing

# 1. When the festival (category 1, one draw) happens. Only the event
#    changes, so the run without events runs once for all values.
# sweep_values = [10, 20, 30, 40, 50, 60, 70]
sweep_label = "event day"
sweep(day) = (; event_config = with_date(event_config, 1, day))

# 2. The average transmission probability (std kept at 0.1, as in
#    1_UserConfig.jl). This changes the whole epidemic, so every value gets
#    its own run without events: twice the runs of sweep 1.
# sweep_values = [0.2, 0.25, 0.3, 0.35, 0.4]
# sweep_label = "mean transmission"
# sweep(m) = (; transmission_dist = Gamma(gamma_params(m, 0.1)...))

# Base population: a .jld2 file holding a table under "data", with at least
# id, age, sex (1 or 2) and the GEMS setting ids (household, and optionally
# office, schoolclass, municipality), like the GEMS population files.
population_file = joinpath(@__DIR__, "Datastorage", "people_Saalekreis.jld2")
# Its matching GEMS settings file (adds e.g. household coordinates, region
# codes and schools/workplaces), or `nothing` if the population has none.
settings_file = nothing # ! Test
#settings_file = joinpath(@__DIR__, "Datastorage", "settings_Saalekreis.jld2") # True settings, slower

scenario = (;
    transmission_dist,
    configfile,
    population = population_file,
    settings_file,
    n_simulations,
    gems_seed,
    run_validation,
    vary_per_replicate,
    run_info_files = [joinpath(@__DIR__, f) for f in ("1_UserConfig.jl", "2_Interface.jl", "config_concert_covid.toml")],
    run_info = (; population_file))

results = if !isnothing(sweep_values)
    run_sweep(event_config; scenario..., sweep, sweep_values, sweep_label)
elseif compare_without_events
    run_comparison(event_config; scenario...)
else
    run_scenario(event_config; scenario...)
end


println("End 2_Interface")
