# ===========================================================================
# scenario.jl
#
# Running a scenario from start to finish (run_scenario): the setup's random
# streams, the population, the GEMS runs (as one batch, or one by one for
# per-replicate setups, branching, picked histories and staged categories),
# and then the analysis, saving and validation.
# ===========================================================================

## === 7. Running a scenario ===

# Where results go unless the caller says otherwise (results_dir in
# run_scenario, run_comparison and run_sweep): the Results/ folder of the
# project this code is in, i.e. next to src/. run.jl and test/test.jl pass
# their own folder explicitly, so they don't depend on where src/ is.
default_results_dir() = normpath(joinpath(@__DIR__, "..", "Results"))

# The setup's random streams, all derived from the setup seed (the TOML's
# Simulation seed). Each setup step gets its own stream, so changing one
# step doesn't reshuffle the others: e.g. changing an event's size or
# n_draws leaves everyone's transmission_prob unchanged.
#   events     - event dates and sizes
#   traits     - per-person traits (transmission_prob)
#   assign     - who attends what
#   population - a synthetic population (synthetic_population)
# Created fresh on every call, so the same seed always gives the same streams.
function setup_streams(seed::Integer)
    setup_rng = Xoshiro(seed)
    return (events     = Xoshiro(rand(setup_rng, UInt64)),
            traits     = Xoshiro(rand(setup_rng, UInt64)),
            assign     = Xoshiro(rand(setup_rng, UInt64)),
            population = Xoshiro(rand(setup_rng, UInt64)))
end

# The stream replicate r draws its own attendance (kind = :attendance) or
# transmission_prob (kind = :transmission) from, when they vary per
# replicate (vary_per_replicate in run_scenario). Derived from the setup
# seed with draws after setup_streams' four, so those stay unchanged.
function replicate_stream(setup_seed::Integer, kind::Symbol, r::Integer)
    setup_rng = Xoshiro(setup_seed)
    draws = [rand(setup_rng, UInt64) for _ in 1:6]   # 1-4: setup_streams
    base_seed = kind == :transmission ? draws[5] : kind == :attendance ? draws[6] :
                error("Unknown kind $kind: use :attendance or :transmission.")
    return Xoshiro(base_seed + UInt64(r))
end

# Calls f with the prepared population in the form GEMS needs. GEMS only
# loads a settings file when the population is given as a file path, not as
# a table. So with a settings file, the prepared population (with the event
# columns) is written to a temporary file and passed by path, then deleted
# afterwards, even if f fails; reading it back first makes sure no event
# column changed on the way. Without one, the table is passed directly,
# which is faster (no file to write, and none to read per replicate).
# GEMS_COLUMNS (attendance.jl) are the columns GEMS copies onto its own
# individuals, for contacts.jl and transmission.jl to read.
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

# Branching: makes future f of history h continue with its own random
# numbers from day branch_day on. GEMS keeps one random generator per
# thread (sim.rngs); at the end of day branch_day - 1 (GEMS runs the
# stepmod before moving to the next day) they're re-seeded in place, with
# copy!, so anything holding on to one of them sees the new state too.
# Future 1 isn't re-seeded: it's the history's own continuation, the same
# as without branching. The seeds come from the GEMS seed, h, f and the
# thread's slot k, so the same settings always give the same futures
# (Xoshiro hashes its seed, so nearby numbers give unrelated streams).
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

# The base population table from a .jld2 file (the table under "data").
# The population files aren't in the repository (see README, Setup), so a
# missing one stops with a message saying what to do.
function load_population(file::AbstractString)
    isfile(file) || error("Population file not found: $file\n" *
        "The population and settings files are not part of the repository. Put them in the data/ " *
        "folder (or set population_file and settings_file in run.jl to where they are); see README, Setup.")
    return JLD2.load(file)["data"]
end

# A synthetic population generated by GEMS, as a table: `synthetic` holds
# the settings of GEMS's generator (n, avg_household_size, avg_office_size,
# avg_school_size); the random numbers come from the setup seed's own
# stream (setup_streams), so the same seed always gives the same population.
# GEMS calls these populations non-representative: they're for trying the
# code out and for tests, not for realistic results.
synthetic_population(synthetic, setup_seed::Integer) =
    dataframe(Population(; synthetic..., rng = setup_streams(setup_seed).population))

# The run_scenario settings `kwargs` with a population file replaced by its
# table, loaded once. Used by run_comparison and run_sweep, which call
# run_scenario several times: without it, every run would read the file
# again (about 10 s for Saalekreis). Same results either way, since
# run_scenario always works on its own copy of the table.
function population_loaded_once(kwargs)
    settings = values(kwargs)   # a NamedTuple
    haskey(settings, :population) && settings.population isa AbstractString || return settings
    return merge(settings, (; population = load_population(settings.population)))
end

# Runs one scenario from start to finish: samples the events, prepares the
# population and assigns attendees, runs the GEMS batch, then analyses,
# saves and (optionally) validates the results. run.jl and test/test.jl
# both call this, so the steps exist only here.
#
#   event_config      - the event categories (config/user_config.jl)
#   transmission_dist - each person's transmission_prob is drawn from it
#   configfile        - the GEMS TOML
#   population        - the base population: a .jld2 file holding a table
#                       under "data", or the table itself, or nothing (the
#                       default) for a synthetic population generated by
#                       GEMS (see synthetic). Never changed: every call
#                       works on its own fresh copy.
#   synthetic         - with population = nothing, the settings of GEMS's
#                       synthetic population generator: n (size),
#                       avg_household_size, avg_office_size, avg_school_size;
#                       default (; n = 100_000). It's made from the setup
#                       seed's own stream, so the same seed gives the same
#                       population.
#   settings_file     - the matching GEMS settings file, or nothing (always
#                       nothing for a synthetic population)
#   n_simulations     - number of replicates
#   setup_seed        - seed of the setup streams; defaults to the TOML's
#                       Simulation seed
#   gems_seed         - GEMS's seed: nothing for a new random one, or a number
#   run_validation    - whether to run validation.jl's checks
#   events_on         - false: everything is set up as usual, but no event
#                       takes place (see daily_global_setting!, contacts.jl)
#   vary_per_replicate - what each replicate draws for itself instead of
#                       sharing: [] (nothing, the default), [:attendance],
#                       [:transmission] (transmission_prob) or both. The
#                       events themselves (dates, sizes) are always shared.
#   histories         - which replicates to run, by number: nothing (all of
#                       1:n_simulations) or e.g. [7], which runs replicate 7
#                       exactly as in a run of 7 or more replicates with the
#                       same seeds and settings
#   branch_day, n_futures - branching: every history runs n_futures times,
#                       identical up to day branch_day - 1 and with its own
#                       random numbers from branch_day on (see branch!).
#                       branch_day = nothing: no branching
#   label             - the GEMS batch label
#   results_dir       - the folder results go into; by default the
#                       project's Results/ (see default_results_dir)
#   results_name      - folder name inside results_dir; by default
#                       run_<n>sims_<date and time>
#   run_info_files    - files copied into run_info/ (see save_run_info)
#   run_info          - extra settings written to run_info.toml
#
# Returns a NamedTuple with everything the run produced: events, people,
# attendees, core_groups (replicate 1's, when they vary), setups (one per
# replicate when something varies, otherwise one shared by all: each with
# attendees, core_groups and transmission_prob), bd (the GEMS BatchData),
# event_results, aggregated, stage_passes (per run: (person id, day) =>
# stage, for infections at later stages of staged categories), run_folder,
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

    # histories (which replicates) and branching (futures per history)
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
    # staged categories need each run's record of which stage infections
    # happened in, so their runs go one by one too
    staged = any(c -> c isa StagedCategory, event_config.categories)

    transmission_func = SettingRate()

    # A GEMS batch of n_runs replicates on the prepared population `pop`
    # (a table, or the file it was written to), whose event columns are in
    # the table `people`.
    make_batch(pop, n_runs, people) = Batch(
        n_runs = n_runs,
        configfile = configfile,
        population = pop,
        settingsfile = settings_file,
        ind_extension = GEMS_COLUMNS,
        transmission_function = transmission_func,
        label = label,
        # GlobalSetting holds only each day's event attendees (see contacts.jl);
        # total event attendances, checked again inside GEMS by daily_global_setting!
        setup = let expected_attendances = sum(length, people.event_dates)
            sim -> daily_global_setting!(sim, expected_attendances; events_on)
        end
    )

    # Picked here rather than by GEMS so the value is known and can be saved
    # with the results (GEMS doesn't report a seed it chose itself).
    gems_seed_used = isnothing(gems_seed) ? rand(RandomDevice(), 0:typemax(Int64)) : Int64(gems_seed)
    Threads.nthreads() == 1 &&
        println("\nRunning on 1 thread. Starting Julia with more, e.g. `julia -t 8 ...`, makes the " *
                "simulation about twice as fast (see README, \"Speed\").")
    println("\nGEMS seed: $gems_seed_used, $(Threads.nthreads()) thread(s). " *
            "To repeat this run, set gems_seed = $gems_seed_used and use the same number of threads.")

    # --- population, attendance and GEMS ---
    # prepare_population and assign_events! add columns to the table, so each
    # call works on a fresh copy and the caller's table stays as it was
    isnothing(population) && !isnothing(settings_file) &&
        error("A settings file belongs to a population file: with population = nothing (a synthetic " *
              "population) settings_file must be nothing too.")
    base = isnothing(population) ? synthetic_population(synthetic, setup_seed) :
           population isa AbstractString ? load_population(population) : copy(population)

    # per run: (person id, day) => stage, for infections in the extra
    # GlobalSetting passes of staged categories (see daily_global_setting!)
    stage_passes = Dict{Tuple{Int32,Int16},Int}[]
    if isempty(vary) && !labelled && !staged
        # one setup, shared by every replicate: one GEMS batch
        people, age_groups, sex_levels = prepare_population(event_config, streams.traits, base, transmission_dist)
        attendees, core_groups = assign_events!(people, events, event_config, streams.assign)
        bd = with_population_for_gems(people, settings_file) do pop
            BatchData(make_batch(pop, n_simulations, people); rd_style = "EssentialResultData", seed = gems_seed_used)
        end
        setups = [(; attendees, core_groups, transmission_prob = people.transmission_prob)]
        append!(stage_passes, [Dict{Tuple{Int32,Int16},Int}() for _ in 1:n_simulations])
    else
        # The runs go one by one, the way GEMS's own batch runs them (the same
        # seed for replicate h as in a normal run), collected into one
        # BatchData. Needed when every replicate draws its own attendance
        # and/or transmission_prob (GEMS then gets each its own population),
        # for branching, and for picking histories. Only the slim part of each
        # setup is kept, to save memory.
        sim_seeds = let rng = Xoshiro(gems_seed_used)
            [GEMS.gems_rand(rng, 0:typemax(Int64)) for _ in 1:maximum(history_list)]
        end
        # the setup of history h: shared by all, or its own (vary_per_replicate)
        function make_setup(h)
            fresh = setup_streams(setup_seed)   # streams from their start, as for replicate-shared setups
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
            # the first history's full table is kept for saving and validation (ids are the same in all)
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
    # replicate r's setup (the same one for all replicates when nothing varies)
    setup_of(r) = setups[length(setups) == 1 ? 1 : r]

    # --- per-event results ---
    event_results = Vector{Dict{String,EventCounts}}()
    for (r, rd) in enumerate(runs(bd))
        push!(event_results, analyze_event_population(infections(rd), people, events, setup_of(r).attendees;
                                                      passes = stage_passes[r]))
    end
    aggregated = aggregate_event_results(event_results)

    # --- saving ---
    # always the project's Results folder, wherever Julia was started from
    name = something(results_name,
        "run_$(length(run_labels))sims_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))")
    run_folder = normpath(joinpath(results_dir, name))
    mkpath(run_folder)

    # everything needed to rerun this run (see README, "Rerunning a run")
    save_run_info(run_folder, run_info_files;
        n_simulations, run_validation, settings_file, events_on,
        population = isnothing(population) ? "synthetic $(synthetic)" :
                     population isa AbstractString ? population : "a table passed in",
        vary_per_replicate = string.(vary), branch_day, n_futures, histories, run_info...,
        setup_seed, gems_seed = gems_seed_used)
    # the per-event results before they're summarised
    save_raw_results(run_folder, events, event_results, setups, people;
                     labels = labelled ? run_labels : nothing)

    println("\nSimulation complete.")
    run_analysis(; aggregated, event_results, events, event_config, bd, run_folder,
                 setup_seed, gems_seed_used, last_day)

    # validation output goes to validation.txt, then is shown on screen too
    if run_validation
        run_validation_checks(; people, events, event_config, attendees,
                              core_groups, aggregated, event_results, bd, run_folder, setups)
    end

    return (; events, people, attendees, core_groups, setups, bd, event_results, aggregated,
              stage_passes, run_folder, setup_seed, gems_seed_used, last_day)
end


