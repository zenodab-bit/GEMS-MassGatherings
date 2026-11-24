# ===========================================================================
# run.jl
#
# Runs a simulation: julia -t 8 run.jl (from this folder). The events are
# in config/user_config.jl; how to run them is set below.
# ===========================================================================

ENV["GKSwstype"] = "100"   # headless plotting; must be set before Plots loads

import Pkg
Pkg.activate(@__DIR__)

using MassGatherings
using Distributions   # for transmission_dist in config/user_config.jl

include("config/user_config.jl")

# === Run ===
n_simulations = 3
run_validation = true
# What each replicate draws for itself instead of sharing (the events
# themselves are always shared): [], [:attendance], [:transmission] or both.
vary_per_replicate = []

# Branching: each replicate runs n_futures times, identical up to day
# branch_day - 1. branch_day = nothing: no branching.
branch_day = nothing
n_futures = 1
# Replicates to run, e.g. [7] to rerun replicate 7 of an earlier run (see
# docs/reproducibility.md); nothing for all.
histories = nothing
# true: also run without the events (same seeds) and compare. Twice as long.
compare_without_events = false

# Sweep: the comparison once per value in sweep_values (nothing: no sweep).
# `sweep` returns what changes for a value. Uncomment one of the two:
sweep_values = nothing

# 1. The festival's day (category 1, one draw).
# sweep_values = [10, 20, 30, 40, 50, 60, 70]
sweep_label = "event day"
sweep(day) = (; event_config = with_date(event_config, 1, day))

# 2. The mean transmission probability (std 0.1). Changes the whole
#    epidemic, so it needs twice the runs of sweep 1.
# sweep_values = [0.2, 0.25, 0.3, 0.35, 0.4]
# sweep_label = "mean transmission"
# sweep(m) = (; transmission_dist = Gamma(gamma_params(m, 0.1)...))

# Population: nothing for a synthetic one (settings in `synthetic`), or a
# .jld2 file, e.g. joinpath(@__DIR__, "data", "people_X.jld2"); see README.
# Keep your own in a copy of this file named my_*.jl (git ignores it).
population_file = nothing
synthetic = (; n = 100_000)
# The population file's GEMS settings file, or nothing (slower).
settings_file = nothing

scenario = (;
    transmission_dist,
    configfile,
    population = population_file,
    synthetic,
    settings_file,
    n_simulations,
    gems_seed,
    run_validation,
    vary_per_replicate,
    branch_day,
    n_futures,
    histories,
    # copied into run_info/
    run_info_files = [joinpath(@__DIR__, "config", "user_config.jl"), @__FILE__, configfile],
    results_dir = joinpath(@__DIR__, "Results"),
    run_info = (; population_file))

results = if !isnothing(sweep_values)
    run_sweep(event_config; scenario..., sweep, sweep_values, sweep_label)
elseif compare_without_events
    run_comparison(event_config; scenario...)
else
    run_scenario(event_config; scenario...)
end


println("End run")
