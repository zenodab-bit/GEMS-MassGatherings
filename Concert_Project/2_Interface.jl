# ===========================================================================
# 2_Interface.jl
#
# Run this file to run a simulation. Loads packages, includes every other
# file, samples events and population from 1_UserConfig.jl's event_config,
# runs the GEMS batch, then triggers analysis and (optionally) validation.
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

events = sample_events(event_config, rng_events)
# sorted for readable output; assign_events! orders dates itself
sort!(events, by = e -> (e.date, e.category_id, e.draw_id, e.section_id))

println("\n=== Sampled Events ===")
for e in events
    println("Event $(e.id) ($(e.name)) — date: $(e.date), n: $(e.n), contacts: $(e.mean_contacts)")
end

sim_length = config["Simulation"]["StopCriterion"]["parameters"]["limit"]
last_day = sim_length - 1   # GEMS simulates days 0 to limit - 1

for e in events
    if e.date > last_day
        @warn "Event $(e.id) on day $(e.date) is after the last simulated day ($last_day; " *
              "GEMS simulates days 0 to StopCriterion limit - 1). Results will be inaccurate."
    end
end

people, age_groups, sex_levels = prepare_population(event_config, rng_traits, JLD2.load(population_file)["data"], transmission_dist)
attendees, core_groups = assign_events!(people, events, event_config, rng_assign)

transmission_func = SettingRate(
    event_dates = Set{Int16}(Int16(e.date) for e in events)
)

# columns GEMS copies onto its own individuals, for 4_Contacts.jl and
# 5_Transmission.jl to read (the list is defined once, in 3_Population.jl)
ind_extension = GEMS_COLUMNS

# GEMS only loads a settings file when the population is given as a file
# path, not as a table. So with a settings file, the prepared population
# (with the event columns) is written to a temporary file and passed by
# path, then deleted after the runs; reading it back first makes sure no
# event column changed on the way. Without one, the table is passed
# directly, which is faster (no file to write, and none to read per replicate).
if isnothing(settings_file)
    population_for_gems = people
else
    population_for_gems = tempname() * ".jld2"
    JLD2.save(population_for_gems, "data", people)
    let saved = JLD2.load(population_for_gems, "data")
        for col in ind_extension
            isequal(saved[!, col], people[!, col]) ||
                error("Column $col changed when writing the population to $population_for_gems.")
        end
    end
end

# total event attendances, checked again inside GEMS by daily_global_setting!
expected_attendances = sum(length, people.event_dates)

b = Batch(
    n_runs = n_simulations,
    configfile = joinpath(@__DIR__, "config_concert_covid.toml"),
    population = population_for_gems,
    settingsfile = settings_file,
    ind_extension = ind_extension,
    transmission_function = transmission_func,
    label = "MultiEvent simulation",
    # GlobalSetting holds only each day's event attendees (see 4_Contacts.jl)
    setup = sim -> daily_global_setting!(sim, expected_attendances)
)
# Picked here rather than by GEMS so the value is known and can be saved
# with the results (GEMS doesn't report a seed it chose itself).
isnothing(gems_seed) || gems_seed isa Integer ||
    error("gems_seed in 1_UserConfig.jl must be nothing or an integer, got $(repr(gems_seed)).")
gems_seed_used = isnothing(gems_seed) ? rand(RandomDevice(), 0:typemax(Int64)) : Int64(gems_seed)
println("\nGEMS seed: $gems_seed_used, $(Threads.nthreads()) thread(s). " *
        "To repeat this run, set gems_seed = $gems_seed_used in 1_UserConfig.jl and use the same number of threads.")

bd = try
    BatchData(b; rd_style = "EssentialResultData", seed = gems_seed_used)
finally
    # delete the temporary population file (if one was written), even if the run fails
    population_for_gems isa String && rm(population_for_gems; force=true)
end

event_results = Vector{Dict{String,EventCounts}}()
for rd in runs(bd)
    inf_log = infections(rd)
    result = analyze_event_population(inf_log, people, events, attendees)
    push!(event_results, result)
end

aggregated = aggregate_event_results(event_results)

# always the project's Results folder, wherever Julia was started from
run_folder = normpath(joinpath(@__DIR__, "..", "Results",
    "run_$(n_simulations)sims_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))"))
mkpath(run_folder)

# everything needed to rerun this run (see README, "Rerunning a run")
save_run_info(run_folder,
    [joinpath(@__DIR__, f) for f in ("1_UserConfig.jl", "2_Interface.jl", "config_concert_covid.toml")];
    n_simulations, run_validation, population_file, settings_file,
    setup_seed = seed, gems_seed = gems_seed_used)
# the per-event results before they're summarised
save_raw_results(run_folder, events, event_results, attendees, core_groups, people)

println("\nSimulation complete.")
run_analysis(; aggregated, event_results, events, event_config, bd, run_folder,
             setup_seed = seed, gems_seed_used, last_day)

# validation output goes to validation.txt, then is shown on screen too
if run_validation
    run_validation_checks(; people, events, event_config, attendees,
                          core_groups, aggregated, bd, run_folder)
end


println("End 2_Interface")