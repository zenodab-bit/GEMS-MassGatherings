# ===========================================================================
# comparison.jl
#
# A scenario with and without its events: run_comparison runs both,
# compare_runs pairs their replicates and writes comparison.txt, .csv, .png.
# ===========================================================================

# Runs a scenario with and without its events, with the same setup and GEMS
# seeds, so replicate r of one pairs with replicate r of the other: each
# pair is identical up to the first event day. Takes run_scenario's
# arguments. Results go to Results/<results_name> (default
# compare_<n>sims_<date and time>): with_events/, without_events/ and the
# comparison files. Validation runs on the run with events only.
# Returns (; with, without, comparison, folder).
function run_comparison(event_config::EventConfig;
                        n_simulations::Integer,
                        gems_seed = nothing,
                        run_validation::Bool = true,
                        results_dir::AbstractString = default_results_dir(),
                        results_name = nothing,
                        kwargs...)

    isnothing(gems_seed) || gems_seed isa Integer ||
        error("gems_seed must be nothing or an integer, got $(repr(gems_seed)).")
    # picked once, so both runs use it
    gems_seed_used = isnothing(gems_seed) ? rand(RandomDevice(), 0:typemax(Int64)) : Int64(gems_seed)

    name = something(results_name,
        "compare_$(n_simulations)sims_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))")
    folder = normpath(joinpath(results_dir, name))

    kwargs = population_loaded_once(kwargs)

    println("\n######## Run 1 of 2: with events ########")
    with = run_scenario(event_config; kwargs..., n_simulations, gems_seed = gems_seed_used, results_dir,
                        run_validation, results_name = joinpath(name, "with_events"))

    println("\n######## Run 2 of 2: without events ########")
    without = run_scenario(event_config; kwargs..., n_simulations, gems_seed = gems_seed_used, results_dir,
                           run_validation = false, events_on = false,
                           results_name = joinpath(name, "without_events"))

    comparison = compare_runs(with, without, folder)
    return (; with, without, comparison, folder)
end


# === Comparison with and without events ===

# Daily new infections on days 0 to last_day, without the seed infections
# (the same in both runs).
function daily_infections(inf_log::DataFrame, last_day::Integer)
    daily = zeros(Int, last_day + 1)
    for (t, s) in zip(inf_log.tick, inf_log.setting_type)
        s == '?' && continue
        0 <= t <= last_day && (daily[t + 1] += 1)
    end
    return daily
end

# The infections before day d, sorted, to compare two runs directly.
infections_before(inf_log::DataFrame, d::Integer) =
    sort!([(t, a, b, s) for (t, a, b, s) in zip(inf_log.tick, inf_log.id_a, inf_log.id_b, inf_log.setting_type) if t < d])

# Compares the runs pair by pair; saves comparison.csv (one row per pair),
# comparison.txt (also shown on screen) and comparison.png into folder, and
# returns the per-pair table.
function compare_runs(with, without, folder)
    events, last_day = with.events, with.last_day
    first_day = isempty(events) ? nothing : minimum(e.date for e in events)
    # per replicate (or one set for all), everyone assigned to an event
    attendee_sets = [Set{Int32}(with.people.id[i] for idxs in values(s.attendees) for i in idxs)
                     for s in with.setups]
    attendee_ids(r) = attendee_sets[length(attendee_sets) == 1 ? 1 : r]

    runs_with, runs_without = runs(with.bd), runs(without.bd)
    length(runs_with) == length(runs_without) ||
        error("The two runs have different numbers of replicates.")

    rows = []
    daily_with, daily_without = Vector{Vector{Int}}(), Vector{Vector{Int}}()
    for (r, (rd_w, rd_o)) in enumerate(zip(runs_with, runs_without))
        log_w, log_o = infections(rd_w), infections(rd_o)
        d_w, d_o = daily_infections(log_w, last_day), daily_infections(log_o, last_day)
        push!(daily_with, d_w)
        push!(daily_without, d_o)

        att_w = count(in(attendee_ids(r)), log_w.id_b)
        att_o = count(in(attendee_ids(r)), log_o.id_b)
        push!(rows, (
            replicate = r,
            identical_before_first_event = isnothing(first_day) ||
                infections_before(log_w, first_day) == infections_before(log_o, first_day),
            total_with = nrow(log_w), total_without = nrow(log_o),
            infected_at_events = count(==('g'), log_w.setting_type),
            attendees_with = att_w, attendees_without = att_o,
            others_with = nrow(log_w) - att_w, others_without = nrow(log_o) - att_o,
            peak_with = maximum(d_w), peak_without = maximum(d_o),
            peak_day_with = argmax(d_w) - 1, peak_day_without = argmax(d_o) - 1,
        ))
    end
    df = DataFrame(identity.(rows))   # identity. gives the rows one concrete type
    mkpath(folder)
    CSV.write(joinpath(folder, "comparison.csv"), df)

    write_comparison_summary(df, with, first_day, joinpath(folder, "comparison.txt"))
    plot_comparison(daily_with, daily_without, events, joinpath(folder, "comparison.png"))
    return df
end

# Per metric: each run's mean and the paired difference (with minus
# without), with a 95% interval for the mean difference (mean ± 1.96 std /
# √pairs). An interval including 0 can't tell the runs apart from chance.
function write_comparison_summary(df, with, first_day, filename)
    n = nrow(df)
    n_identical = count(df.identical_before_first_event)

    function row(io, label, w, o)
        d = w .- o
        sd = n > 1 ? std(d) : NaN
        half = 1.96 * sd / sqrt(n)
        ci = n > 1 ? "$(fmt(mean(d) - half)) to $(fmt(mean(d) + half))" : "n/a"
        println(io, rpad(label, 34), " | ", rpad(fmt(mean(w)), 9), " | ", rpad(fmt(mean(o)), 9), " | ",
                rpad(fmt(mean(d)), 9), " | ", rpad(n > 1 ? fmt(sd) : "n/a", 9), " | ",
                rpad(fmt(minimum(d)), 9), " | ", rpad(fmt(maximum(d)), 9), " | ", ci)
    end

    open(filename, "w") do io
        println(io, "=== Comparison: with vs. without events ===")
        println(io, "Replicate pairs: $n")
        println(io, "Setup seed: $(with.setup_seed); GEMS seed: $(with.gems_seed_used) (the same for both runs)")
        if isnothing(first_day)
            println(io, "No events in this scenario: both runs are the same scenario.")
        else
            println(io, "First event day: $first_day")
            println(io, n_identical == n ?
                "✓ All $n pairs are identical before day $first_day, so the differences come from the events." :
                "✗ Only $n_identical of $n pairs are identical before day $first_day: the runs are not paired " *
                "properly, so the differences below are not only the events' effect.")
        end

        println(io, "\nDifference = with events minus without, per pair; 95% interval = for the mean difference.")
        println(io, rpad("Metric", 34), " | ", rpad("With", 9), " | ", rpad("Without", 9), " | ",
                rpad("Diff mean", 9), " | ", rpad("Diff std", 9), " | ", rpad("Diff min", 9), " | ",
                rpad("Diff max", 9), " | ", "95% interval")
        println(io, "-"^125)
        row(io, "Total infections", df.total_with, df.total_without)
        row(io, "Infections among attendees", df.attendees_with, df.attendees_without)
        row(io, "Infections among everyone else", df.others_with, df.others_without)
        row(io, "Peak daily infections", df.peak_with, df.peak_without)
        row(io, "Peak day", df.peak_day_with, df.peak_day_without)

        println(io, "\nInfected at events (direct count, run with events only): $(fmt(mean(df.infected_at_events))) on average")
        # per pair: extra infections overall per infection at the events
        ratios = [(r.total_with - r.total_without) / r.infected_at_events
                  for r in eachrow(df) if r.infected_at_events > 0]
        if isempty(ratios)
            println(io, "Extra infections per infection at events: n/a (no infections at events)")
        else
            skipped = n - length(ratios)
            println(io, "Extra infections per infection at events: $(fmt(mean(ratios)))" *
                        (skipped > 0 ? " ($skipped pair(s) with no infections at events skipped)" : ""))
            println(io, "  Above 1: the events cause more infections than happen at them (onward spread).")
            println(io, "  Below 1: some event infections would have happened later anyway.")
        end
        println(io, "\nPer-pair values are in comparison.csv.")
    end

    println("\n", read(filename, String))
    println("Saved: comparison.txt, comparison.csv")
end

# Mean daily new infections with and without events, and the cumulative
# extra infections (mean, range over pairs as a ribbon), event days marked.
function plot_comparison(daily_with, daily_without, events, filename)
    days = 0:(length(daily_with[1]) - 1)
    mat_w, mat_o = hcat(daily_with...), hcat(daily_without...)

    p1 = plot(title="Daily new infections", xlabel="Day", ylabel="New infections",
        legend=:outertopright)
    plot!(p1, days, mean(mat_w, dims=2)[:], label="With events", color=PLOT_COLORS.vermillion, linewidth=2)
    plot!(p1, days, mean(mat_o, dims=2)[:], label="Without events", color=PLOT_COLORS.blue, linewidth=2)
    add_event_vlines!(p1, events)

    extra = cumsum(mat_w .- mat_o, dims=1)
    avg = mean(extra, dims=2)[:]
    lo, hi = minimum(extra, dims=2)[:], maximum(extra, dims=2)[:]
    p2 = plot(title="Extra infections so far (with minus without)", xlabel="Day",
        ylabel="Cumulative difference", legend=:outertopright)
    plot!(p2, days, avg, ribbon=(avg .- lo, hi .- avg), fillalpha=0.2,
        label="Mean (range over pairs)", color=PLOT_COLORS.green, linewidth=2)
    hline!(p2, [0], color=:gray, linestyle=:dot, label="")
    add_event_vlines!(p2, events)

    p = plot(p1, p2, layout=(2, 1), size=(1000, 800), dpi=300; PLOT_MARGINS...)
    savefig(p, filename)
    println("Saved: comparison.png")
end

