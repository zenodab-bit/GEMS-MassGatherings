# ===========================================================================
# 6_Analysis.jl
#
# Auto-run at the end of 2_Interface.jl (after a successful simulation) to
# write text metrics and generate every plot. Not meant to be included on
# its own — the calls at the bottom of this file (print_metrics,
# plot_epidemic_overview, etc.) execute immediately using variables
# (aggregated, events, bd, run_folder, event_results, n_simulations)
# assumed to already exist in the calling scope.
#
# fmt, print_table_header, print_metric_row, and add_event_vlines!, used
# throughout this file, are defined in 0_Helpers.jl (included first).
# ===========================================================================


# === Text metrics ===

# Writes per-event compartment metrics and global epidemic stats to
# multievent_analysis.txt in run_folder.
function print_metrics(aggregated, events, bd, run_folder)
    filename = "$run_folder/multievent_analysis.txt"

    open(filename, "w") do io
        println(io, "=== MultiEvent Analysis ===")
        println(io, "Simulations: $n_simulations")
        println(io, "Events: $(length(events))")

        # --- Global epidemic metrics ---
        println(io, "\n--- Global Epidemic Metrics ---")
        println(io, "Total infections: $(fmt(total_infections(bd)["mean"]))")
        println(io, "Attack rate:      $(fmt(attack_rate(bd)["mean"]))")
        println(io, "R0:               $(fmt(r0(bd)["mean"]))")

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
            print_metric_row(io, "Dead", metrics[:dead])
            print_metric_row(io, "Same day other", metrics[:same_day_other])
            print_metric_row(io, "Infected at event", metrics[:infected_at_event])

            infection_rate = metrics[:infected_at_event].mean / event.n * 100
            println(io, "\nInfection rate: $(fmt(infection_rate))%")
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
        'h' => :orange,
        'c' => :red,
        'o' => :purple,
        'g' => :blue,
        'm' => :brown,
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
    fields = [:susceptible, :infectious, :exposed, :recovered, :dead]
    colors = [:green, :red, :blue, :gray, :black]

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
function plot_infected_boxplot(event_results, events, run_folder)
    event_ids = [e.id for e in events]
    values_per_event = [[r[e.id].infected_at_event for r in event_results] for e in events]

    p = plot(title="Infected at Event", xlabel="Event", ylabel="Infected", dpi=300,
        legend=:outertopright)

    if length(event_results) > 1
        for (i, vals) in enumerate(values_per_event)
            boxplot!(p, [event_ids[i]], vals, label=event_ids[i], legend=true)
        end
    else
        single_vals = [vals[1] for vals in values_per_event]
        bar!(p, event_ids, single_vals, label=false)
    end

    savefig(p, "$run_folder/infected_boxplot.png")
    println("Saved: infected_boxplot")
end

# Like the old sitting/standing boxplot: one panel per (category, section),
# one box per draw at its date. rate=true plots % of attendees instead of
# counts, since n differs between draws. Saved to infected_by_section.png.
function plot_infected_by_section(event_results, events, event_config, run_folder; rate=true)
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
            v = Float64[r[e.id].infected_at_event for r in event_results]
            rate && (v .*= 100 / e.n)
            append!(vals, v)
            append!(days, fill(e.date, length(v)))
        end

        # boxplot needs >1 replicate to mean anything; single run falls back to points
        plotfn = length(event_results) > 1 ? boxplot : scatter
        push!(panels, plotfn(days, vals, title=title, xlabel="Event day",
            ylabel=rate ? "Infected (% of attendees)" : "Infected",
            xticks=unique(days), legend=false))
    end

    p = plot(panels..., layout=(length(panels), 1), size=(800, 350 * length(panels)), dpi=300)
    savefig(p, "$run_folder/infected_by_section.png")
    println("Saved: infected_by_section")
end

# === Auto-generate ===
print_metrics(aggregated, events, bd, run_folder)
plot_epidemic_overview(bd, events, run_folder)
plot_cases_by_setting(bd, events, run_folder)
plot_event_seir(aggregated, events, run_folder)
plot_infected_boxplot(event_results, events, run_folder)
plot_infected_by_section(event_results, events, event_config, run_folder)

println("\nAnalysis complete. Results saved to $run_folder")
println("Analysis")