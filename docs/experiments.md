# Experiments: comparisons, sweeps, per-replicate setups, branching

All set in `run.jl` (or your `my_*.jl` copy) and combinable.

## Comparing with and without events

Infections *at* events aren't the infections events *cause*: some would have happened later anyway, and onward infections elsewhere aren't counted. With `compare_without_events = true` the scenario runs twice: with events, and with the same setup but nobody going. Both use the same seeds, so replicate *r* of one pairs with replicate *r* of the other; each pair is identical up to the first event day, so later differences come from the events. `comparison.txt` checks this.

Results go to `Results/compare_<n>sims_<date and time>/`:

| File | Content |
|---|---|
| `with_events/`, `without_events/` | The two runs (validation on `with_events/` only) |
| `comparison.txt` | Per metric (infections overall, among attendees and others, peak height and day): each run's mean and the paired difference, with a 95% interval (including 0: run more replicates). Also the extra infections per infection at events (above 1: onward spread; below 1: some would have happened anyway). |
| `comparison.csv` | The values per pair |
| `comparison.png` | Daily infections with and without events, and the cumulative difference |

A few events can have a small effect next to the variation between runs, so comparisons may need more replicates. Pairs drift apart after the first event, so totals over long runs are noisiest. With a high attack rate, events mostly change *when* people are infected. A comparison takes twice as long.

## Sweeps

A sweep runs the comparison once per value of one setting:

- `sweep_values`: the values, or `nothing` (default);
- `sweep_label`: their name in the output;
- `sweep`: a function returning the settings that change for a value.

`run.jl` has two ready-made sweeps; uncomment the one you want:

```julia
# 1. when the festival (category 1, one draw) happens
sweep_values = [10, 20, 30, 40, 50, 60, 70]
sweep_label = "event day"
sweep(day) = (; event_config = with_date(event_config, 1, day))

# 2. the average transmission probability, with the std kept at 0.1
sweep_values = [0.2, 0.25, 0.3, 0.35, 0.4]
sweep_label = "mean transmission"
sweep(m) = (; transmission_dist = Gamma(gamma_params(m, 0.1)...))
```

`with_date` needs `n_draws = 1`; for recurring categories, `with_date_range(event_config, id, (first, last))` moves the window. Both work for staged categories.

All values use the same seeds. If `sweep` only changes `event_config`, the run without events is the same for every value and runs once: 7 values cost 8 runs instead of 14.

Results go to `Results/sweep_<n>sims_<date and time>/`: a `value_<value>/` folder per value (its runs and comparison), the shared `without_events/`, and `sweep.txt`, `.csv` and `.png` with one row or point per value.

## Varying the setup per replicate

By default all replicates share one setup (events, `transmission_prob`, attendees), so they show the epidemic's variation *given that draw*. If the draw happens to be unusual, every replicate inherits it, and more replicates don't change that. `vary_per_replicate` lets each replicate draw parts for itself:

| `vary_per_replicate` | Each replicate draws its own | Its spread then also shows |
|---|---|---|
| `[]` (default) | nothing | (only the epidemic's own randomness) |
| `[:attendance]` | who attends what | how much it matters *who* goes |
| `[:transmission]` | everyone's `transmission_prob` | how much it matters how infectious individuals happen to be |
| `[:attendance, :transmission]` | both | both |

Event dates and sizes stay shared. Each replicate's draw comes from its own stream, so the same seeds give the same draws. `attendance.jld2` then holds one attendee list per replicate, and validation checks replicate 1's assignment. Replicates run one by one, a few seconds slower each.

## Branching

```julia
branch_day = 44   # the futures differ from this day on (nothing: no branching)
n_futures = 20    # futures per history
```

Every replicate (a **history**) runs `n_futures` times, identical up to day `branch_day - 1`; then GEMS's random generators are re-seeded so each **future** goes its own way. Future 1 is the history itself. The seeds come from `gems_seed`, the history and the future, so the same settings give the same futures. Use `n_simulations = 1` for futures of one situation.

**Picking histories**: `histories = [7]` runs only replicate 7, exactly as in the full run with the same `gems_seed` and settings (see [Rerunning a run](reproducibility.md#rerunning-a-run)). Several work too, e.g. `[3, 7]`.

Every future counts as a replicate; `event_results.csv` gets `history` and `future` columns. Comparisons and sweeps pair the same history and future.
