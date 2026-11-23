# ===========================================================================
# utils.jl
#
# Small helpers used throughout: contact-number and probability distributions
# (negative binomial, gamma, beta), age groups, the plot colours, and table
# formatting for the text output.
# ===========================================================================

## === 1. Numeric distribution helpers ===

# (r, p) params for NegativeBinomial. Requires variance > mean.
# use sample_n_contacts, not this directly, for inputs that might violate that.
function negbin_params(mean, std)
    variance = std^2
    p = mean / variance
    r = mean^2 / (variance - mean)
    return r, p
end

# Draws a contact count from (mean, std); falls back to Poisson(mean) when
# variance <= mean, since NegativeBinomial isn't valid there.
function sample_n_contacts(rng, mean, std)
    mean <= 0 && return 0
    variance = std^2
    if variance <= mean
        return rand(rng, Poisson(mean))
    end
    r, p = negbin_params(mean, std)
    return rand(rng, NegativeBinomial(r, p))
end

# (α, β) params for Gamma. Valid for any mean, std > 0, but Gamma is
# unbounded above — clamp to [0,1] if the result is used as a probability
# (see prepare_population in attendance.jl).
function gamma_params(mean, std)
    variance = std^2
    β = variance / mean
    α = mean / β
    return α, β
end

# (α, β) params for Beta, which always stays within [0, 1]. Needs
# 0 < mean < 1 and std^2 < mean * (1 - mean); errors otherwise.
function beta_params(mean, std)
    0 < mean < 1 || error("beta_params: mean ($mean) must be between 0 and 1.")
    variance = std^2
    variance < mean * (1 - mean) ||
        error("beta_params: std ($std) too large for mean $mean; a Beta needs std^2 < mean * (1 - mean) = $(mean * (1 - mean)).")
    k = mean * (1 - mean) / variance - 1
    return mean * k, (1 - mean) * k
end


## === 2. Age-group helpers ===

# Classifies a raw age into a group index (1..length(age_boundaries)+1) via
# binary search. age_boundaries must be sorted ascending.
age_group_idx(age, age_boundaries) = searchsortedfirst(age_boundaries, age)

# Formats a group index as a label, e.g. index 1 of [45,65] → "<=45".
function age_group_label_from_idx(idx::Int, age_boundaries)
    idx > length(age_boundaries) && return ">$(age_boundaries[end])"
    return "<=$(age_boundaries[idx])"
end

# Classify and format a raw age in one call.
age_group_label(age, age_boundaries) = age_group_label_from_idx(age_group_idx(age, age_boundaries), age_boundaries)



## === 6. Output formatting & plotting helpers ===

# Colours for every plot: the Okabe-Ito palette, chosen to stay
# distinguishable for the common forms of colour blindness. Change a colour
# here to change it everywhere. analysis.jl makes PLOT_PALETTE the default
# palette, so plots that don't set colours themselves (e.g. GEMS's epidemic
# overview) use it in this order. Black is left out of the palette, since it
# hides the median lines in boxplots and clashes with the event-day lines.
const PLOT_COLORS = (
    blue       = colorant"#0072B2",
    orange     = colorant"#E69F00",
    green      = colorant"#009E73",
    vermillion = colorant"#D55E00",
    purple     = colorant"#CC79A7",
    skyblue    = colorant"#56B4E9",
    yellow     = colorant"#F0E442",   # low contrast on white: last in the palette
    black      = colorant"#000000",
)
const PLOT_PALETTE = [PLOT_COLORS.blue, PLOT_COLORS.orange, PLOT_COLORS.green, PLOT_COLORS.vermillion,
                      PLOT_COLORS.purple, PLOT_COLORS.skyblue, PLOT_COLORS.yellow]

# Rounds a number for consistent display in printed tables.
function fmt(x; digits=2)
    round(x, digits=digits)
end

# Column-header row for the per-event metrics table (multievent_analysis.txt).
function print_table_header(io)
    println(io, rpad("Metric", 20), " | ",
        rpad("Mean", 9), " | ", rpad("Std", 9), " | ",
        rpad("Min", 9), " | ", rpad("P25", 9), " | ",
        rpad("Median", 9), " | ", rpad("P75", 9), " | ", "Max")
    println(io, "-"^100)
end

# One row of the per-event metrics table.
function print_metric_row(io, label, m)
    println(io, rpad(label, 20), " | ",
        rpad(fmt(m.mean), 9), " | ", rpad(fmt(m.std), 9), " | ",
        rpad(fmt(m.min), 9), " | ", rpad(fmt(m.p25), 9), " | ",
        rpad(fmt(m.median), 9), " | ", rpad(fmt(m.p75), 9), " | ",
        fmt(m.max))
end

# One dashed vline per unique event date, labeled with every event.id
# sharing that date (avoids duplicate overlapping lines/legend entries).
function add_event_vlines!(p, events)
    seen_dates = Set{Int}()
    for event in events
        event.date in seen_dates && continue
        push!(seen_dates, event.date)

        same_date_ids = [e.id for e in events if e.date == event.date]
        label = length(same_date_ids) > 1 ?
                "Events $(join(same_date_ids, ", ")) (day $(event.date))" :
                "Event $(event.id) (day $(event.date))"

        vline!(p, [event.date],
            linestyle=:dash,
            color=:black,
            label=label,
            linewidth=1.5)
    end
end

println("End utils")
