# ===========================================================================
# utils.jl
#
# Small helpers: contact-number and probability distributions, age groups,
# plot colours and table formatting.
# ===========================================================================

## === Distributions ===

# (r, p) for NegativeBinomial; needs variance > mean (use sample_n_contacts).
function negbin_params(mean, std)
    variance = std^2
    p = mean / variance
    r = mean^2 / (variance - mean)
    return r, p
end

# A contact count with the given mean and std: negative binomial, or
# Poisson(mean) when variance <= mean, where the negative binomial isn't valid.
function sample_n_contacts(rng, mean, std)
    mean <= 0 && return 0
    variance = std^2
    if variance <= mean
        return rand(rng, Poisson(mean))
    end
    r, p = negbin_params(mean, std)
    return rand(rng, NegativeBinomial(r, p))
end

# (α, β) for Gamma, for any mean, std > 0. Unbounded above, so clamp it
# when used as a probability (as prepare_population does).
function gamma_params(mean, std)
    variance = std^2
    β = variance / mean
    α = mean / β
    return α, β
end

# (α, β) for Beta, which stays within [0, 1]. Needs 0 < mean < 1 and
# std^2 < mean * (1 - mean).
function beta_params(mean, std)
    0 < mean < 1 || error("beta_params: mean ($mean) must be between 0 and 1.")
    variance = std^2
    variance < mean * (1 - mean) ||
        error("beta_params: std ($std) too large for mean $mean; a Beta needs std^2 < mean * (1 - mean) = $(mean * (1 - mean)).")
    k = mean * (1 - mean) / variance - 1
    return mean * k, (1 - mean) * k
end


## === Age groups ===

# Group index of an age (1..length(age_boundaries)+1); boundaries ascending.
age_group_idx(age, age_boundaries) = searchsortedfirst(age_boundaries, age)

# Label of a group index, e.g. 1 of [45,65] → "<=45".
function age_group_label_from_idx(idx::Int, age_boundaries)
    idx > length(age_boundaries) && return ">$(age_boundaries[end])"
    return "<=$(age_boundaries[idx])"
end

age_group_label(age, age_boundaries) = age_group_label_from_idx(age_group_idx(age, age_boundaries), age_boundaries)



## === Output formatting and plots ===

# Plot colours: the Okabe-Ito palette, readable with common colour
# blindness. run_analysis makes PLOT_PALETTE the default palette (GEMS's
# plots use it too). Black isn't in it: it hides boxplot medians and
# clashes with the event-day lines.
const PLOT_COLORS = (
    blue       = colorant"#0072B2",
    orange     = colorant"#E69F00",
    green      = colorant"#009E73",
    vermillion = colorant"#D55E00",
    purple     = colorant"#CC79A7",
    skyblue    = colorant"#56B4E9",
    yellow     = colorant"#F0E442",   # low contrast on white: last
    black      = colorant"#000000",
)
const PLOT_PALETTE = [PLOT_COLORS.blue, PLOT_COLORS.orange, PLOT_COLORS.green, PLOT_COLORS.vermillion,
                      PLOT_COLORS.purple, PLOT_COLORS.skyblue, PLOT_COLORS.yellow]

# Rounds a number for printed tables.
function fmt(x; digits=2)
    round(x, digits=digits)
end

# Header row of the per-event table in multievent_analysis.txt.
function print_table_header(io)
    println(io, rpad("Metric", 20), " | ",
        rpad("Mean", 9), " | ", rpad("Std", 9), " | ",
        rpad("Min", 9), " | ", rpad("P25", 9), " | ",
        rpad("Median", 9), " | ", rpad("P75", 9), " | ", "Max")
    println(io, "-"^100)
end

# One row of that table.
function print_metric_row(io, label, m)
    println(io, rpad(label, 20), " | ",
        rpad(fmt(m.mean), 9), " | ", rpad(fmt(m.std), 9), " | ",
        rpad(fmt(m.min), 9), " | ", rpad(fmt(m.p25), 9), " | ",
        rpad(fmt(m.median), 9), " | ", rpad(fmt(m.p75), 9), " | ",
        fmt(m.max))
end

# One dashed line per event date, labelled with every event on it.
function add_event_vlines!(p, events)
    seen_dates = Set{Int}()
    for event in events
        event.date in seen_dates && continue
        push!(seen_dates, event.date)

        # a count, not every id: a long label squeezes the plot (and makes GKS complain)
        n_same_date = count(e -> e.date == event.date, events)
        label = n_same_date > 1 ?
                "Day $(event.date): $n_same_date events" :
                "Event $(event.id) (day $(event.date))"

        vline!(p, [event.date],
            linestyle=:dash,
            color=:black,
            label=label,
            linewidth=1.5)
    end
end

