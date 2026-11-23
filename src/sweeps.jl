# ===========================================================================
# sweeps.jl
#
# Sweeps: the comparison with and without events for a list of values of
# one setting (run_sweep), the date helpers to build them (with_date,
# with_date_range), and the summary over all values (sweep.txt, sweep.csv,
# sweep.png).
# ===========================================================================

## === 8. Sweeps ===

# A copy of event_config with category `category_id`'s date_range set to
# `date_range`; everything else unchanged. For sweeps over when events happen.
function with_date_range(event_config::EventConfig, category_id::Integer, date_range::Tuple{Int,Int})
    i = findfirst(c -> c.id == category_id, event_config.categories)
    i === nothing && error("No category with id $category_id in event_config.")
    categories = copy(event_config.categories)
    c = categories[i]
    categories[i] = typeof(c)(c; date_range)   # Category or StagedCategory, the rest unchanged
    return EventConfig(event_config; categories)
end

# The same, with every draw of the category on one day. Only for categories
# with n_draws = 1: draws of one category can't share a day (use
# with_date_range to move the window of a recurring category instead).
function with_date(event_config::EventConfig, category_id::Integer, day::Integer)
    i = findfirst(c -> c.id == category_id, event_config.categories)
    i === nothing && error("No category with id $category_id in event_config.")
    n = event_config.categories[i].n_draws
    n <= 1 || error("Category $category_id has n_draws = $n, and its draws can't all be on day $day. " *
                    "Use with_date_range to move its date window instead.")
    return with_date_range(event_config, category_id, (Int(day), Int(day)))
end

# Runs the comparison with and without events (run_comparison) once per
# value in sweep_values, to see how the events' effect depends on one
# setting. `sweep` says what changes: sweep(value) returns the run_scenario
# settings to use for that value, as a NamedTuple, e.g.
#   sweep(day) = (; event_config = with_date(event_config, 1, day))
#   sweep(m)   = (; transmission_dist = Gamma(gamma_params(m, 0.1)...))
# Everything not returned stays as passed in (kwargs, the same as for
# run_scenario). sweep_label names the setting in the output.
#
# Every value uses the same setup seed and the same GEMS seed, so values
# differ only in the swept setting, not in luck.
#
# Shortcut: if sweep only ever returns event_config, the run without events
# is the same for every value (nobody goes to any event, and the
# population, transmission_prob and seeds don't change), so it runs once
# and every value is compared against it. Otherwise each value gets its own.
#
# Results go to Results/<results_name>, by default
# sweep_<n>sims_<date and time>: one folder per value (value_<v>/, holding
# with_events/, its comparison files and, without the shortcut, its own
# without_events/), the shared without_events/ (with the shortcut), and the
# summary over all values: sweep.txt, sweep.csv, sweep.png.
#
# Returns (; summary, folder): the per-value summary table and the folder.
function run_sweep(event_config::EventConfig;
                   sweep::Function,
                   sweep_values::AbstractVector,
                   sweep_label::AbstractString = "value",
                   n_simulations::Integer,
                   gems_seed = nothing,
                   run_validation::Bool = true,
                   results_name = nothing,
                   kwargs...)

    isempty(sweep_values) && error("sweep_values is empty.")
    allunique(sweep_values) || error("sweep_values has repeated values: $(sweep_values).")
    isnothing(gems_seed) || gems_seed isa Integer ||
        error("gems_seed must be nothing or an integer, got $(repr(gems_seed)).")
    # picked once here, so every run uses the same one
    gems_seed_used = isnothing(gems_seed) ? rand(RandomDevice(), 0:typemax(Int64)) : Int64(gems_seed)

    kwargs = population_loaded_once(kwargs)
    overrides = [sweep(v) for v in sweep_values]
    for (v, o) in zip(sweep_values, overrides)
        o isa NamedTuple ||
            error("sweep($v) returned $(typeof(o)); it must return settings like (; event_config = ...).")
    end
    events_only = all(o -> keys(o) == (:event_config,), overrides)

    name = something(results_name,
        "sweep_$(n_simulations)sims_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))")
    folder = normpath(joinpath(@__DIR__, "..", "Results", name))
    n_runs = length(sweep_values) + (events_only ? 1 : length(sweep_values))
    println("\n=== Sweep: $(length(sweep_values)) values of $sweep_label, $n_runs runs of " *
            "$n_simulations replicate(s) ($(events_only ? "one run without events for all values" :
                                              "one run without events per value")) ===")

    # settings shared by every run; for each value, its own settings win
    shared = (; kwargs..., n_simulations, gems_seed = gems_seed_used)
    settings_for(o) = merge(shared, Base.structdiff(o, NamedTuple{(:event_config,)}))
    config_for(o) = get(o, :event_config, event_config)

    shared_without = nothing
    if events_only
        println("\n######## Without events (shared by every value) ########")
        shared_without = run_scenario(config_for(overrides[1]); settings_for(overrides[1])...,
            run_validation = false, events_on = false,
            results_name = joinpath(name, "without_events"))
    end

    rows = []
    for (v, o) in zip(sweep_values, overrides)
        value_name = joinpath(name, "value_$v")
        println("\n######## $sweep_label = $v: with events ########")
        with = run_scenario(config_for(o); settings_for(o)..., run_validation,
                            results_name = joinpath(value_name, "with_events"))
        without = if events_only
            shared_without
        else
            println("\n######## $sweep_label = $v: without events ########")
            run_scenario(config_for(o); settings_for(o)..., run_validation = false, events_on = false,
                         results_name = joinpath(value_name, "without_events"))
        end
        df = compare_runs(with, without, normpath(joinpath(@__DIR__, "..", "Results", value_name)))
        push!(rows, summarise_comparison(v, df))
    end

    summary = DataFrame(identity.(rows))   # identity. gives the rows one concrete type
    write_sweep_summary(summary, sweep_label, events_only, gems_seed_used, folder)
    return (; summary, folder)
end


# === Sweeps ===
# Used by run_sweep (above): one comparison per value, summarised
# into one row per value.

# Mean of x and its 95% interval (mean ± 1.96 std / √n); NaN bounds for n = 1.
function mean_interval(x)
    m = mean(x)
    length(x) > 1 || return (mean=m, lo=NaN, hi=NaN)
    half = 1.96 * std(x) / sqrt(length(x))
    return (mean=m, lo=m - half, hi=m + half)
end

# One value's comparison table (compare_runs) as one summary row.
function summarise_comparison(value, df)
    extra = mean_interval(df.total_with .- df.total_without)
    ratios = [(r.total_with - r.total_without) / r.infected_at_events
              for r in eachrow(df) if r.infected_at_events > 0]
    return (
        value = value,
        pairs = nrow(df),
        pairs_identical_before_events = count(df.identical_before_first_event),
        infected_at_events = mean(df.infected_at_events),
        extra_infections = extra.mean, extra_lo = extra.lo, extra_hi = extra.hi,
        extra_attendees = mean(df.attendees_with .- df.attendees_without),
        extra_others = mean(df.others_with .- df.others_without),
        extra_per_event_infection = isempty(ratios) ? NaN : mean(ratios),
        peak_change = mean(df.peak_with .- df.peak_without),
        peak_day_change = mean(df.peak_day_with .- df.peak_day_without),
        peak_day_without = mean(df.peak_day_without),
    )
end

# Saves sweep.csv (one row per value), sweep.txt (the same as a table, also
# shown on screen) and sweep.png into folder.
function write_sweep_summary(summary, sweep_label, events_only, gems_seed_used, folder)
    mkpath(folder)
    CSV.write(joinpath(folder, "sweep.csv"), summary)

    filename = joinpath(folder, "sweep.txt")
    open(filename, "w") do io
        println(io, "=== Sweep over $sweep_label ===")
        println(io, "Values: $(join(summary.value, ", ")); $(summary.pairs[1]) replicate pair(s) each")
        println(io, "GEMS seed: $gems_seed_used (the same for every run)")
        println(io, events_only ?
            "One run without events, shared by every value (only the events change between values)." :
            "One run without events per value (the swept setting changes more than the events).")
        all(summary.pairs_identical_before_events .== summary.pairs) ||
            println(io, "✗ Some pairs are not identical before the first event day: see the values' " *
                        "comparison.txt. Their differences are not only the events' effect.")

        println(io, "\nPer value: extra infections = with events minus without (mean over pairs, " *
                    "95% interval for the mean); the other columns are means over pairs.")
        cols = [("Infected at events", :infected_at_events), ("Extra infections", :extra_infections),
                ("Extra: attendees", :extra_attendees), ("Extra: others", :extra_others),
                ("Extra per event inf.", :extra_per_event_infection), ("Peak change", :peak_change),
                ("Peak day change", :peak_day_change)]
        println(io, rpad(sweep_label, 12), " | ", join([rpad(c[1], 20) for c in cols], " | "), " | 95% interval (extra)")
        println(io, "-"^(15 + 23 * length(cols) + 20))
        for r in eachrow(summary)
            interval = isnan(r.extra_lo) ? "n/a" : "$(fmt(r.extra_lo)) to $(fmt(r.extra_hi))"
            println(io, rpad(string(r.value), 12), " | ",
                    join([rpad(string(fmt(r[c[2]])), 20) for c in cols], " | "), " | ", interval)
        end
        println(io, "\nExtra per event infection: above 1, the events cause more infections than happen " *
                    "at them (onward spread); below 1, some event infections would have happened later anyway.")
        println(io, "Peak day without events, for reference: $(join(unique(fmt.(summary.peak_day_without)), ", ")).")
        println(io, "Each value's own results and comparison are in value_<value>/.")
    end
    println("\n", read(filename, String))

    plot_sweep(summary, sweep_label, joinpath(folder, "sweep.png"))
    println("Saved: sweep.txt, sweep.csv, sweep.png")
end

# Two panels against the swept value: extra infections (with the 95%
# interval) next to the infections counted at events, and the extra
# infections per event infection. Saved to filename.
function plot_sweep(summary, sweep_label, filename)
    x = summary.value
    p1 = plot(title="Effect of the events", xlabel=sweep_label, ylabel="Infections",
        legend=:outertopright)
    has_interval = !any(isnan, summary.extra_lo)
    plot!(p1, x, summary.extra_infections, label="Extra infections (with minus without)",
        color=PLOT_COLORS.green, linewidth=2, marker=:circle,
        yerror=has_interval ? (summary.extra_infections .- summary.extra_lo,
                               summary.extra_hi .- summary.extra_infections) : nothing)
    plot!(p1, x, summary.infected_at_events, label="Infected at events (direct count)",
        color=PLOT_COLORS.vermillion, linewidth=2, marker=:square)
    hline!(p1, [0], color=:gray, linestyle=:dot, label="")

    p2 = plot(title="Extra infections per infection at events", xlabel=sweep_label,
        ylabel="Ratio", legend=false)
    plot!(p2, x, summary.extra_per_event_infection, color=PLOT_COLORS.blue, linewidth=2, marker=:circle)
    hline!(p2, [1], color=:gray, linestyle=:dot)

    p = plot(p1, p2, layout=(2, 1), size=(1000, 800), dpi=300)
    savefig(p, filename)
end

println("End sweeps")
