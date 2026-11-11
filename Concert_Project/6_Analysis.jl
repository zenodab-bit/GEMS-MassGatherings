# ===========================================================================
# 6_Analysis.jl
#
# Text metrics and plots for a finished run. Including this file only
# defines functions; run_analysis(...) at the bottom runs them all.
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
        println(io, "\n--- Global Epidemic Metrics ---")
        println(io, "Total infections: $(fmt(total_infections(bd)["mean"]))")
        println(io, "Attack rate:      $(fmt(attack_rate(bd)["mean"]))")
        println(io, "R0:               $(fmt(r0(bd)["mean"]))")

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

println("End Analysis")