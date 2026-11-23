# ===========================================================================
# 6_Analysis.jl
#
# Text metrics and plots for a finished run. Including this file only
# defines functions; run_analysis(...) runs them all. At the bottom,
# compare_runs compares a run with events against the same run without
# (see run_comparison in 0_Helpers.jl), and write_sweep_summary sums up a
# sweep (see run_sweep in 0_Helpers.jl).
#
# fmt, print_table_header, print_metric_row, add_event_vlines!,
# susceptible_attack_rates, offspring_dispersion, PLOT_COLORS and
# PLOT_PALETTE, used throughout this file, are defined in 0_Helpers.jl
# (included first).
# ===========================================================================


# === Text metrics ===

# Writes per-event compartment metrics and global epidemic stats to
# multievent_analysis.txt in run_folder.
function print_metrics(aggregated, event_results, events, bd, run_folder, setup_seed, gems_seed_used, last_day)
    filename = "$run_folder/multievent_analysis.txt"

    open(filename, "w") do io
        println(io, "=== MultiEvent Analysis ===")
        println(io, "Simulations: $(length(event_results))")
        println(io, "Events: $(length(events))")
        println(io, "Setup seed (TOML): $setup_seed")
        println(io, "GEMS seed: $gems_seed_used ($(Threads.nthreads()) thread(s); " *
                    "set gems_seed to this value in the config you ran with and use the same thread count to repeat this run)")

        # --- Global epidemic metrics ---
        # GEMS.-qualified, so a variable of the same name elsewhere (e.g. r0) can't hide them
        println(io, "\n--- Global Epidemic Metrics ---")
        println(io, "Total infections: $(fmt(GEMS.total_infections(bd)["mean"]))")
        println(io, "Attack rate:      $(fmt(GEMS.attack_rate(bd)["mean"]))")
        println(io, "R0:               $(fmt(GEMS.r0(bd)["mean"]))")

        # dispersion of onward infections, one value per run
        disp = [offspring_dispersion(infections(rd), last_day) for rd in runs(bd)]
        ks = [d.k for d in disp if !isnan(d.k)]
        k_str = isempty(ks) ? "n/a (too few infections)" :
                all(isinf, ks) ? "Inf (no overdispersion)" :
                "$(fmt(median(ks))) (median over runs; min $(fmt(minimum(ks))), max $(fmt(maximum(ks))))"
        cutoffs = [d.cutoff for d in disp]
        cutoff_str = minimum(cutoffs) == maximum(cutoffs) ? "$(cutoffs[1])" :
                     "$(minimum(cutoffs))-$(maximum(cutoffs))"
        println(io, "Offspring k:      $k_str")
        println(io, "                  lower k = more superspreading, Inf = no overdispersion. People infected in the " *
                    "last $cutoff_str days are left out: 95% of infected people stopped being able to infect within " *
                    "that time in these runs (until recovery, death or hospital admission)")

        # --- Per event metrics ---
        for event in events
            println(io, "\n\n##############################")
            println(io, "### Event $(event.id) — Day $(event.date) — n=$(event.n)")
            println(io, "##############################")

            print_table_header(io)
            metrics = aggregated[event.id]
            print_metric_row(io, "Susceptible", metrics[:susceptible])
            print_metric_row(io, "Infectious", metrics[:infectious])
            print_metric_row(io, "Exposed", metrics[:exposed])
            print_metric_row(io, "Recovered", metrics[:recovered])
            print_metric_row(io, "Absent", metrics[:absent])   # dead, severe or hospitalised that day
            print_metric_row(io, "Same day other", metrics[:same_day_other])
            print_metric_row(io, "Infected at event", metrics[:infected_at_event])

            # per run, among attendees who were actually there (assigned
            # minus absent); runs where nobody was there are skipped
            shares = [100 * r[event.id].infected_at_event / (event.n - r[event.id].absent)
                      for r in event_results if event.n > r[event.id].absent]
            share_str = isempty(shares) ? "n/a (no attendees present)" : "$(fmt(mean(shares)))%"
            println(io, "\nShare of present attendees infected:    $share_str")

            # mean of per-run rates; runs with no susceptible attendees are skipped
            rates = susceptible_attack_rates(event_results, event.id)
            rate_str = isempty(rates) ? "n/a (no susceptible attendees)" : "$(fmt(mean(rates)))%"
            skipped = length(event_results) - length(rates)
            skipped > 0 && (rate_str *= " ($skipped run(s) with no susceptible attendees skipped)")
            println(io, "Attack rate among susceptible attendees: $rate_str")
        end
    end

    println("Metrics saved to $filename")
end


# === Plots ===

# Three-panel plot (tick cases, cumulative cases, effective R) with event
# dates marked, saved to epidemic_overview.png.
function plot_epidemic_overview(bd, events, run_folder)
    p1 = gemsplot(bd, type=:TickCases)
    add_event_vlines!(p1, events)

    p2 = gemsplot(bd, type=:CumulativeCases)
    add_event_vlines!(p2, events)

    p3 = gemsplot(bd, type=:EffectiveReproduction)
    add_event_vlines!(p3, events)

    p_overview = plot(p1, p2, p3, layout=(3, 1), size=(800, 900), dpi=300,
        legend=:outertopright)
    savefig(p_overview, "$run_folder/epidemic_overview.png")
    println("Saved: epidemic_overview")
end


# Daily case counts by setting type (mean with min/max ribbon across
# replicates), saved to cases_by_setting.png. Excludes the 'Initial'
# setting_type ('?') — a single-point seeding spike at tick 0, not a trend.
function plot_cases_by_setting(bd, events, run_folder)
    setting_colors = Dict(
        'h' => PLOT_COLORS.blue,
        'c' => PLOT_COLORS.green,
        'o' => PLOT_COLORS.purple,
        'g' => PLOT_COLORS.vermillion,   # events: the strongest colour
        'm' => PLOT_COLORS.skyblue,
    )
    setting_labels = Dict(
        'h' => "Household",
        'c' => "SchoolClass",
        'o' => "Office",
        'g' => "GlobalSetting",
        'm' => "Municipality",
    )

    all_setting_cases = Dict{Char,Vector{Vector{Float64}}}()

    # same length for every replicate, so they can be stacked even when one
    # run's last infection comes earlier than another's
    n_ticks = maximum(maximum(infections(rd).tick) for rd in runs(bd))
    for rd in runs(bd)
        inf_log = infections(rd)
        tick_setting = combine(groupby(inf_log, [:tick, :setting_type]), nrow => :daily_cases)

        for s in keys(setting_labels)
            rows = filter(r -> r.setting_type == s, tick_setting)
            series = zeros(Float64, n_ticks + 1)
            for row in eachrow(rows)
                series[row.tick+1] = Float64(row.daily_cases)
            end
            if !haskey(all_setting_cases, s)
                all_setting_cases[s] = Vector{Vector{Float64}}()
            end
            push!(all_setting_cases[s], series)
        end
    end

    p = plot(title="Cases by Setting", xlabel="Tick", ylabel="Daily Cases", dpi=300,
        legend=:outerright, size=(1400, 600))
    for (s, series_list) in all_setting_cases
        isempty(series_list) && continue
        mat = hcat(series_list...)
        avg = mean(mat, dims=2)[:]
        lo = minimum(mat, dims=2)[:]
        hi = maximum(mat, dims=2)[:]
        any(avg .> 0) || continue
        plot!(p, 0:(length(avg)-1), avg,
            ribbon=(avg .- lo, hi .- avg),
            fillalpha=0.2,
            label=setting_labels[s],
            color=setting_colors[s],
            linewidth=2)
    end
    add_event_vlines!(p, events)
    savefig(p, "$run_folder/cases_by_setting.png")
    println("Saved: cases_by_setting")
end


# One SEIR-compartment bar chart per event, saved as seir_<event.name>.png.
# Filenames use event.name, not the guaranteed-unique event.id — category
# and section names need to stay distinct for this not to collide.
function plot_event_seir(aggregated, events, run_folder)
    fields = [:susceptible, :infectious, :exposed, :recovered, :absent]   # absent: dead, severe or hospitalised
    colors = [PLOT_COLORS.skyblue, PLOT_COLORS.vermillion, PLOT_COLORS.orange,
              PLOT_COLORS.green, PLOT_COLORS.black]

    for event in events
        metrics = aggregated[event.id]
        values = [metrics[f].mean for f in fields]

        p = bar(string.(fields), values,
            color=colors,
            title="SEIR State — $(event.name) Day $(event.date)",
            xlabel="Compartment",
            ylabel="Count",
            legend=false,
            dpi=300)
        savefig(p, "$run_folder/seir_$(event.name).png")
        println("Saved: seir_$(event.name)")
    end
end


# Boxplot of infected-at-event across replicates, one box per event, saved
# to infected_boxplot.png. Falls back to a plain bar chart with a single
# value per event when there's only one replicate (n_simulations == 1),
# since a boxplot needs a distribution to be meaningful.
# Boxes are coloured by category (one legend entry per category); each box
# is labelled with its event id on the x axis.
function plot_infected_boxplot(event_results, events, event_config, run_folder)
    # one palette colour per category, in config order
    cat_color = Dict(c.id => PLOT_PALETTE[mod1(i, length(PLOT_PALETTE))]
                     for (i, c) in enumerate(event_config.categories))
    cat_name = Dict(c.id => (isempty(c.name) ? "category $(c.id)" : c.name)
                    for c in event_config.categories)

    p = plot(title="Infected at Event", xlabel="Event", ylabel="Infected", dpi=300,
        legend=:outertopright)

    labelled = Set{Int}()   # categories already in the legend
    for e in events
        vals = [r[e.id].infected_at_event for r in event_results]
        label = e.category_id in labelled ? "" : cat_name[e.category_id]
        push!(labelled, e.category_id)
        if length(event_results) > 1
            boxplot!(p, [e.id], vals, color=cat_color[e.category_id], label=label)
        else
            bar!(p, [e.id], vals, color=cat_color[e.category_id], label=label)
        end
    end

    savefig(p, "$run_folder/infected_boxplot.png")
    println("Saved: infected_boxplot")
end

# Like the old sitting/standing boxplot: one panel per (category, section),
# one box per draw at its date. rate=true plots the attack rate among
# susceptible attendees (per run, as in multievent_analysis.txt) instead of
# counts, since n and the number of susceptibles differ between draws; runs
# with no susceptible attendees are skipped. Saved to infected_by_section.png.
function plot_infected_by_section(event_results, events, event_config, run_folder; rate=true)
    isempty(events) && return
    panels = []
    for (cat_id, sec_id) in unique((e.category_id, e.section_id) for e in events)
        category = event_config.categories[findfirst(c -> c.id == cat_id, event_config.categories)]
        section = category.sections[findfirst(s -> s.id == sec_id, category.sections)]
        cat_label = isempty(category.name) ? "category $cat_id" : category.name
        sec_label = !isempty(section.name) ? section.name :
                    length(category.sections) > 1 ? "sec$sec_id" : ""
        title = isempty(sec_label) ? cat_label : "$cat_label ($sec_label)"

        days = Int[]
        vals = Float64[]
        for e in filter(e -> e.category_id == cat_id && e.section_id == sec_id, events)
            v = rate ? susceptible_attack_rates(event_results, e.id) :
                       Float64[r[e.id].infected_at_event for r in event_results]
            append!(vals, v)
            append!(days, fill(e.date, length(v)))
        end

        if isempty(vals)
            push!(panels, plot(title="$title: no susceptible attendees", legend=false, axis=false, grid=false))
            continue
        end

        # boxplot needs >1 replicate to mean anything; single run falls back to points
        plotfn = length(event_results) > 1 ? boxplot : scatter
        push!(panels, plotfn(days, vals, title=title, xlabel="Event day",
            ylabel=rate ? "Infected (% of susceptible attendees)" : "Infected",
            xticks=unique(days), legend=false))
    end

    p = plot(panels..., layout=(length(panels), 1), size=(800, 350 * length(panels)), dpi=300)
    savefig(p, "$run_folder/infected_by_section.png")
    println("Saved: infected_by_section")
end

# === Run everything ===

# Writes the text metrics and every plot for a finished run into run_folder.
function run_analysis(; aggregated, event_results, events, event_config, bd,
                        run_folder, setup_seed, gems_seed_used, last_day)
    # colour-blind-friendly palette for every plot that doesn't set its own
    # colours, including GEMS's (see PLOT_COLORS in 0_Helpers.jl)
    default(palette = PLOT_PALETTE)

    print_metrics(aggregated, event_results, events, bd, run_folder, setup_seed, gems_seed_used, last_day)
    plot_epidemic_overview(bd, events, run_folder)
    plot_cases_by_setting(bd, events, run_folder)
    plot_event_seir(aggregated, events, run_folder)
    plot_infected_boxplot(event_results, events, event_config, run_folder)
    plot_infected_by_section(event_results, events, event_config, run_folder)

    println("\nAnalysis complete. Results saved to $run_folder")
end


# === Comparison with and without events ===
# Used by run_comparison (0_Helpers.jl). Replicate r of the run with events
# is paired with replicate r of the run without: same seeds, so the pair is
# identical up to the first event day and every difference after it comes
# from the events.

# Daily new infections on days 0 to last_day. Initial (seed) infections are
# left out: they're the same in both runs.
function daily_infections(inf_log::DataFrame, last_day::Integer)
    daily = zeros(Int, last_day + 1)
    for (t, s) in zip(inf_log.tick, inf_log.setting_type)
        s == '?' && continue
        0 <= t <= last_day && (daily[t + 1] += 1)
    end
    return daily
end

# The infections before day d, sorted, so two runs can be compared directly.
infections_before(inf_log::DataFrame, d::Integer) =
    sort!([(t, a, b, s) for (t, a, b, s) in zip(inf_log.tick, inf_log.id_a, inf_log.id_b, inf_log.setting_type) if t < d])

# Compares the two runs pair by pair. Saves comparison.csv (one row per
# pair), comparison.txt (the summary, also shown on screen) and
# comparison.png into folder, and returns the per-pair table.
function compare_runs(with, without, folder)
    events, last_day = with.events, with.last_day
    first_day = isempty(events) ? nothing : minimum(e.date for e in events)
    # per replicate, everyone assigned to at least one event in the run with
    # events (one shared set unless attendance varies per replicate)
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

# The summary table: per metric, the mean of each run and the paired
# difference (with minus without) across pairs. The 95% interval is for the
# mean difference (mean ± 1.96 std / √pairs): if it includes 0, the pairs
# can't tell the difference apart from chance, and more replicates help.
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
        # per pair: extra infections overall per infection counted at the events
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

# Two panels: mean daily new infections with and without events, and the
# cumulative extra infections (with minus without; mean, with the range
# over pairs as a ribbon). Event days are marked. Saved to filename.
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

    p = plot(p1, p2, layout=(2, 1), size=(1000, 800), dpi=300)
    savefig(p, filename)
    println("Saved: comparison.png")
end


# === Sweeps ===
# Used by run_sweep (0_Helpers.jl): one comparison per value, summarised
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

println("End Analysis")