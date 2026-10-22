# ===========================================================================
# 2_Interface.jl
#
# Run this file to run a simulation. Loads packages, includes every other
# file, samples events and population from 1_UserConfig.jl's event_config,
# runs the GEMS batch, then triggers analysis and (optionally) validation.
#
# The two settings below (n_simulations, run_validation) are the main
# things to adjust here; everything else you'd want to change lives in
# 1_UserConfig.jl.
# ===========================================================================

ENV["GKSwstype"] = "100"   # headless plotting backend; must be set before Plots loads

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GEMS, Parameters, DataFrames, Distributions, CSV,
    CategoricalArrays, JLD2, Random, StatsBase, Plots, StatsPlots, Dates, TOML
using Random: Xoshiro, shuffle

# === Includes ===
include("0_Helpers.jl")
include("1_UserConfig.jl")
include("3_Population.jl")
include("4_Contacts.jl")
include("5_Transmission.jl")

# === Run ===
n_simulations = 3
run_validation = true

events = sample_events(event_config, rng)
# sorted for readable output; assign_events! orders dates itself
sort!(events, by = e -> (e.date, e.category_id, e.draw_id, e.section_id))

println("\n=== Sampled Events ===")
for e in events
    println("Event $(e.id) ($(e.name)) — date: $(e.date), n: $(e.n), contacts: $(e.mean_contacts)")
end

sim_length = config["Simulation"]["StopCriterion"]["parameters"]["limit"]

for e in events
    if e.date > sim_length
        @warn "Event $(e.id) on day $(e.date) is beyond simulation length ($sim_length). Results will be inaccurate."
    end
end

people, age_groups, sex_levels = prepare_population(event_config, rng)
attendees, core_groups = assign_events!(people, events, event_config, rng)

transmission_func = SettingRate(
    general_rate = general_rate,
    event_dates = Set{Int16}(Int16(e.date) for e in events)
)

b = Batch(
    n_runs = n_simulations,
    configfile = joinpath(@__DIR__, "config_concert_covid.toml"),
    population = people,
    settingsfile = joinpath(@__DIR__, "Datastorage", "settings_Saalekreis.jld2"),
    ind_extension = [:category_ids, :draw_ids, :section_ids, :mean_event_contacts, 
    :std_event_contacts, :event_dates, :transmission_prob, :cross_section_mean_contacts, :cross_section_std_contacts],
    transmission_function = transmission_func,
    label = "MultiEvent simulation"
)
bd = BatchData(b; rd_style = "EssentialResultData")

event_results = Vector{Dict{String,EventCounts}}()
for rd in runs(bd)
    inf_log = infections(rd)
    result = analyze_event_population(inf_log, people, events, attendees)
    push!(event_results, result)
end

aggregated = aggregate_event_results(event_results)

run_folder = "Results/run_$(n_simulations)sims_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM"))"
mkpath(run_folder)

println("\nSimulation complete.")
include("6_Analysis.jl")

if run_validation
    include("7_Validation.jl")
end


println("End 2_Interface")