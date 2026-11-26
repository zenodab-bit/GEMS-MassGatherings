# ===========================================================================
# scenario.jl
#
# run_scenario: one scenario from start to finish. The setup's random
# streams, the population, the GEMS runs (as one batch, or one by one), then
# analysis, saving and validation.
# ===========================================================================

## === Running a scenario ===

# Results/ next to src/, unless the caller passes results_dir (run.jl and
# test/test.jl do).
default_results_dir() = normpath(joinpath(@__DIR__, "..", "Results"))

# One random stream per setup step, from the setup seed, so changing one
# step doesn't reshuffle the others.
function setup_streams(seed::Integer)
    setup_rng = Xoshiro(seed)
    return (events     = Xoshiro(rand(setup_rng, UInt64)),
            traits     = Xoshiro(rand(setup_rng, UInt64)),
            assign     = Xoshiro(rand(setup_rng, UInt64)),
            population = Xoshiro(rand(setup_rng, UInt64)))
end

# Replicate r's own stream for :attendance or :transmission
# (vary_per_replicate), from draws after setup_streams' four.
function replicate_stream(setup_seed::Integer, kind::Symbol, r::Integer)
    setup_rng = Xoshiro(setup_seed)
    draws = [rand(setup_rng, UInt64) for _ in 1:6]   # 1-4: setup_streams
    base_seed = kind == :transmission ? draws[5] : kind == :attendance ? draws[6] :
                error("Unknown kind $kind: use :attendance or :transmission.")
    return Xoshiro(base_seed + UInt64(r))
end

# Calls f with the population as GEMS needs it: the table, or with a
# settings file a temporary file (GEMS reads settings only for those).
function with_population_for_gems(f, people, settings_file)
    isnothing(settings_file) && return f(people)
    file = tempname() * ".jld2"
    try
        JLD2.save(file, "data", people)
        let saved = JLD2.load(file, "data")
            for col in GEMS_COLUMNS
                isequal(saved[!, col], people[!, col]) ||
                    error("Column $col changed when writing the population to $file.")
            end
        end
        return f(file)
    finally
        rm(file; force=true)
    end
end

# Branching: future f > 1 of history h gets its own random numbers from
# branch_day on: GEMS's per-thread generators (sim.rngs) are re-seeded at
# the end of day branch_day - 1. In place, with copy!, so anything else
# holding one of them sees the new state too; replacing them wouldn't do
# that. Future 1 is the history itself. Seeds come from the GEMS seed, h, f
# and the thread slot k, so the same settings give the same futures.
future_seed(gems_seed, h, f, k) =
    reinterpret(UInt64, Int64(gems_seed)) ⊻ (UInt64(h) << 40) ⊻ (UInt64(f) << 20) ⊻ UInt64(k)

function branch!(sim, branch_day::Integer, gems_seed, h::Integer, f::Integer)
    f == 1 && return sim
    reseed!(s) = for (k, g) in enumerate(s.rngs)
        copy!(g, Xoshiro(future_seed(gems_seed, h, f, k)))
    end
    if branch_day == 0
        reseed!(sim)   # before the first day
    else
        previous = sim.stepmod   # keep the stepmod already set (daily_global_setting!)
        sim.stepmod = s -> (previous(s); tick(s) == branch_day - 1 && reseed!(s); nothing)
    end
    return sim
end

# The population table in a .jld2 file (under "data"); population files
# aren't in the repository, so a missing one says what to do.
function load_population(file::AbstractString)
    isfile(file) || error("Population file not found: $file\n" *
        "The population and settings files are not part of the repository. Put them in the data/ " *
        "folder (or set population_file and settings_file in run.jl to where they are); see README, \"The population\".")
    return JLD2.load(file)["data"]
end

# A synthetic population from GEMS's generator (settings in `synthetic`),
# from the setup seed's own stream. Not representative: for trying and tests.
synthetic_population(synthetic, setup_seed::Integer) =
    dataframe(Population(; synthetic..., rng = setup_streams(setup_seed).population))

# The settings `kwargs` with a population file replaced by its table, so
# comparisons and sweeps read it only once (about 10 s for Saalekreis).
function population_loaded_once(kwargs)
    settings = values(kwargs)   # a NamedTuple
    haskey(settings, :population) && settings.population isa AbstractString || return settings
    return merge(settings, (; population = load_population(settings.population)))
end

# Runs one scenario: checks the config, samples the events, prepares the
# population and assigns attendees, runs GEMS, then analyses, saves and
# (optionally) validates. run.jl, test.jl, comparisons and sweeps all call it.
#
#   event_config, transmission_dist, configfile - as in config/user_config.jl
#   population        - a .jld2 file (table under "data"), a table, or
#                       nothing for a synthetic one (settings in
#                       `synthetic`); never changed: each call uses a copy
#   settings_file     - the population file's GEMS settings file, or nothing
#   setup_seed        - default: the TOML's Simulation seed
#   gems_seed         - nothing for a random one (printed and saved)
#   events_on         - false: set up as usual, but no event takes place
#   vary_per_replicate - [], [:attendance], [:transmission] or both: drawn
#                       per replicate instead of shared (event dates and
#                       sizes are always shared)
#   histories         - replicates to run, e.g. [7], exactly as in a full
#                       run; nothing for 1:n_simulations
#   branch_day, n_futures - see branch!; branch_day = nothing: no branching
#   results_dir, results_name - default Results/run_<n>sims_<date and time>
#   run_info_files, run_info - copied to / written in run_info/
#
# Returns a NamedTuple: events, people, attendees and core_groups
# (replicate 1's when they vary), setups (one per replicate when something
# varies, else one), bd (GEMS BatchData), event_results, aggregated,
# stage_passes (per run: (person id, day) => stage), run_folder,
# setup_seed, gems_seed_used and last_day.
function run_scenario(event_config::EventConfig;
                      transmission_dist,
                      configfile::AbstractString,
                      population = nothing,
                      synthetic = (; n = 100_000),
                      settings_file = nothing,
                      n_simulations::Integer,
                      setup_seed::Integer = TOML.parsefile(configfile)["Simulation"]["seed"],
                      gems_seed = nothing,
                      run_validation::Bool = true,
                      events_on::Bool = true,
                      vary_per_replicate = Symbol[],
                      branch_day = nothing,
                      n_futures::Integer = 1,
                      histories = nothing,
                      label::AbstractString = "MultiEvent simulation",
                      results_dir::AbstractString = default_results_dir(),
                      results_name = nothing,
                      run_info_files::Vector{String} = [configfile],
                      run_info = (;))

    isnothing(gems_seed) || gems_seed isa Integer ||
        error("gems_seed must be nothing or an integer, got $(repr(gems_seed)).")
    # here too, so configs built in code (e.g. by a sweep) are checked
    validate_config(event_config; quiet = true)
    streams = setup_streams(setup_seed)

    # --- events ---
    events = sample_events(event_config, streams.events)
    # sorted for readable output; assign_events! orders dates itself
    sort!(events, by = e -> (e.date, e.category_id, e.draw_id, e.section_id))

    println("\n=== Sampled Events ===")
    for e in events
        println("Event $(e.id) ($(e.name)) — date: $(e.date), n: $(e.n), contacts: $(e.mean_contacts)")
    end

    sim_length = TOML.parsefile(configfile)["Simulation"]["StopCriterion"]["parameters"]["limit"]
    last_day = sim_length - 1   # GEMS simulates days 0 to limit - 1

    for e in events
        if e.date > last_day
            @warn "Event $(e.id) on day $(e.date) is after the last simulated day ($last_day; " *
                  "GEMS simulates days 0 to StopCriterion limit - 1). Results will be inaccurate."
        end
    end

    vary = Symbol.(collect(vary_per_replicate))
    all(in((:attendance, :transmission)), vary) && allunique(vary) ||
        error("vary_per_replicate may only contain :attendance and :transmission, got $(repr(vary_per_replicate)).")

    isnothing(histories) ||
        (histories isa AbstractVector && !isempty(histories) && all(h -> h isa Integer && h >= 1, histories) &&
         allunique(histories)) ||
        error("histories must be nothing or a list of different replicate numbers (1, 2, ...), got $(repr(histories)).")
    n_futures >= 1 || error("n_futures must be at least 1, got $n_futures.")
    n_futures == 1 || !isnothing(branch_day) ||
        error("n_futures = $n_futures needs a branch_day: the day from which the futures differ.")
    isnothing(branch_day) || (branch_day isa Integer && 0 <= branch_day <= last_day) ||
        error("branch_day must be nothing or a day from 0 to $last_day (the last simulated day), got $(repr(branch_day)).")
    history_list = isnothing(histories) ? collect(1:n_simulations) : collect(Int, histories)
    # one entry per GEMS run, in run order
    run_labels = [(history = h, future = f) for h in history_list for f in 1:n_futures]
    labelled = !isnothing(branch_day) || !isnothing(histories)
    # staged categories need each run's pass_infections, so they run one by one too
    staged = any(c -> c isa StagedCategory, event_config.categories)

    transmission_func = SettingRate()

    # A GEMS batch of n_runs replicates on `pop` (a table, or its file),
    # whose event columns are in `people`.
    make_batch(pop, n_runs, people) = Batch(
        n_runs = n_runs,
        configfile = configfile,
        population = pop,
        settingsfile = settings_file,
        ind_extension = GEMS_COLUMNS,
        transmission_function = transmission_func,
        label = label,
        # GlobalSetting holds only each day's attendees (contacts.jl)
        setup = let expected_attendances = sum(length, people.event_dates)
            sim -> daily_global_setting!(sim, expected_attendances; events_on)
        end
    )

    # picked here, not by GEMS, so it can be saved with the results
    gems_seed_used = isnothing(gems_seed) ? rand(RandomDevice(), 0:typemax(Int64)) : Int64(gems_seed)
    Threads.nthreads() == 1 &&
        println("\nRunning on 1 thread. Starting Julia with more, e.g. `julia -t 8 ...`, makes the " *
                "simulation about twice as fast (see docs/reproducibility.md).")
    println("\nGEMS seed: $gems_seed_used, $(Threads.nthreads()) thread(s). " *
            "To repeat this run, set gems_seed = $gems_seed_used and use the same number of threads.")

    # --- population, attendance and GEMS ---
    # every call works on a fresh copy; the caller's table stays as it was
    isnothing(population) && !isnothing(settings_file) &&
        error("A settings file belongs to a population file: with population = nothing (a synthetic " *
              "population) settings_file must be nothing too.")
    base = isnothing(population) ? synthetic_population(synthetic, setup_seed) :
           population isa AbstractString ? load_population(population) : copy(population)

    # per run: (person id, day) => stage, for infections in staged categories' extra passes
    stage_passes = Dict{Tuple{Int32,Int16},Int}[]
    if isempty(vary) && !labelled && !staged
        # one setup shared by every replicate: one GEMS batch
        people, age_groups, sex_levels = prepare_population(event_config, streams.traits, base, transmission_dist)
        attendees, core_groups = assign_events!(people, events, event_config, streams.assign)
        bd = with_population_for_gems(people, settings_file) do pop
            BatchData(make_batch(pop, n_simulations, people); rd_style = "EssentialResultData", seed = gems_seed_used)
        end
        setups = [(; attendees, core_groups, transmission_prob = people.transmission_prob)]
        append!(stage_passes, [Dict{Tuple{Int32,Int16},Int}() for _ in 1:n_simulations])
    else
        # one by one, with the seeds GEMS's batch would use, into one
        # BatchData; only the slim part of each setup is kept
        sim_seeds = let rng = Xoshiro(gems_seed_used)
            [GEMS.gems_rand(rng, 0:typemax(Int64)) for _ in 1:maximum(history_list)]
        end
        # history h's setup: shared, or its own (vary_per_replicate)
        function make_setup(h)
            fresh = setup_streams(setup_seed)   # from their start, as for shared setups
            traits_rng = :transmission in vary ? replicate_stream(setup_seed, :transmission, h) : fresh.traits
            assign_rng = :attendance in vary ? replicate_stream(setup_seed, :attendance, h) : fresh.assign
            people_h, _, _ = prepare_population(event_config, traits_rng, copy(base), transmission_dist)
            attendees_h, core_groups_h = assign_events!(people_h, events, event_config, assign_rng)
            return people_h, (; attendees = attendees_h, core_groups = core_groups_h,
                                transmission_prob = people_h.transmission_prob)
        end
        shared = isempty(vary) ? make_setup(1) : nothing
        processor = GEMS.BatchProcessor(; keep_rundata = true, master_seed = gems_seed_used)
        setups = []
        people = nothing
        for h in history_list
            people_h, setup_h = isnothing(shared) ? make_setup(h) : shared
            # the first history's table is kept for saving and validation (same ids in all)
            isnothing(people) && (people = people_h)

            with_population_for_gems(people_h, settings_file) do pop
                for f in 1:n_futures
                    println("\n=== Replicate $h" * (n_futures > 1 ? ", future $f of $n_futures" : "") *
                            (isempty(vary) ? "" : ": its own $(join(string.(vary), " and "))") * " ===")
                    push!(setups, setup_h)
                    b = make_batch(pop, 1, people_h)
                    sim = Simulation(; only(GEMS.simconfigs(b))..., seed = sim_seeds[h])
                    only(GEMS.simsetups(b))(sim)
                    isnothing(branch_day) || branch!(sim, branch_day, gems_seed_used, h, f)
                    run!(sim)
                    push!(stage_passes, copy(settings(sim, GlobalSetting)[1].contact_sampling_method.pass_infections))
                    GEMS.accumulate!(processor, PostProcessor(sim); rd_style = "EssentialResultData")
                end
            end
        end
        bd = BatchData(processor; style = "DefaultBatchData")
        isempty(vary) && (setups = [shared[2]])   # one setup for all runs
        attendees, core_groups = setups[1].attendees, setups[1].core_groups
    end
    # replicate r's setup (the same one for all when nothing varies)
    setup_of(r) = setups[length(setups) == 1 ? 1 : r]

    # --- per-event results ---
    event_results = Vector{Dict{String,EventCounts}}()
    for (r, rd) in enumerate(runs(bd))
        push!(event_results, analyze_event_population(infections(rd), people, events, setup_of(r).attendees;
                                                      passes = stage_passes[r]))
    end
    aggregated = aggregate_event_results(event_results)

    # --- saving ---
    name = something(results_name,
        "run_$(length(run_labels))sims_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))")
    run_folder = normpath(joinpath(results_dir, name))
    mkpath(run_folder)

    # everything needed to rerun this run (docs/reproducibility.md)
    save_run_info(run_folder, run_info_files;
        n_simulations, run_validation, settings_file, events_on,
        population = isnothing(population) ? "synthetic $(synthetic)" :
                     population isa AbstractString ? population : "a table passed in",
        vary_per_replicate = string.(vary), branch_day, n_futures, histories, run_info...,
        setup_seed, gems_seed = gems_seed_used)
    save_raw_results(run_folder, events, event_results, setups, people;
                     labels = labelled ? run_labels : nothing)

    println("\nSimulation complete.")
    run_analysis(; aggregated, event_results, events, event_config, bd, run_folder,
                 setup_seed, gems_seed_used, last_day, events_on)

    # written to validation.txt, then shown on screen
    if run_validation
        run_validation_checks(; people, events, event_config, attendees,
                              core_groups, aggregated, event_results, bd, run_folder, setups)
    end

    return (; events, people, attendees, core_groups, setups, bd, event_results, aggregated,
              stage_passes, run_folder, setup_seed, gems_seed_used, last_day)
end

