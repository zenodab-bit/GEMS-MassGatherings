# ===========================================================================
# 2_Interface.jl
#
# Run this file to run a simulation. Loads packages, includes every other
# file, then calls run_scenario (0_Helpers.jl), which samples events and
# population from 1_UserConfig.jl's event_config, runs the GEMS batch, then
# triggers analysis and (optionally) validation.
#
# The settings below (n_simulations, run_validation and the population
# files) are the main things to adjust here; everything else you'd want to
# change lives in 1_UserConfig.jl.
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

# Base population: a .jld2 file holding a table under "data", with at least
# id, age, sex (1 or 2) and the GEMS setting ids (household, and optionally
# office, schoolclass, municipality), like the GEMS population files.
population_file = joinpath(@__DIR__, "Datastorage", "people_Saalekreis.jld2")
# Its matching GEMS settings file (adds e.g. household coordinates, region
# codes and schools/workplaces), or `nothing` if the population has none.
settings_file = nothing # ! Test
#settings_file = joinpath(@__DIR__, "Datastorage", "settings_Saalekreis.jld2") # True settings, slower

results = run_scenario(event_config;
    transmission_dist,
    configfile,
    population = population_file,
    settings_file,
    n_simulations,
    gems_seed,
    run_validation,
    run_info_files = [joinpath(@__DIR__, f) for f in ("1_UserConfig.jl", "2_Interface.jl", "config_concert_covid.toml")],
    run_info = (; population_file))


println("End 2_Interface")
