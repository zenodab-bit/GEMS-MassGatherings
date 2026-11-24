# Experiments: comparisons, sweeps, per-replicate setups, branching

All of these are set in `run.jl` (or your own `my_*.jl` copy of it) and can be combined.

## Comparing with and without events

The per-event results count who was infected *at* an event. That isn't the same as how many infections the events *cause*: some people infected at an event would have been infected later anyway (at home, at work), and the people they go on to infect at home or at work aren't counted as event infections. To measure the events' real effect, set `compare_without_events = true`. The scenario then runs twice:

1. **with events**, as usual;
2. **without events**: the same events and attendees are set up, but nobody goes, so the events don't take place.

Both runs use the same setup seed and the same GEMS seed (picked once if `gems_seed = nothing`), so replicate 1 of one run pairs with replicate 1 of the other, and so on. Before the first event day, nobody is at an event in either run, so each pair is identical up to that day, and every difference after it comes from the events. `comparison.txt` checks this and says so.

Results go to `Results/compare_<n>sims_<date and time>/`:

| File | Content |
|---|---|
| `with_events/`, `without_events/` | The two runs, each a normal run folder. Validation runs on `with_events/` only, since its checks assume the events took place. In `without_events/`, the per-event results show the would-be attendees on the event day. |
| `comparison.txt` | Per metric (total infections, infections among attendees and among everyone else, peak height and day): the mean of each run, and the difference (with minus without) per pair, summarised as mean, std, min, max and a 95% interval for the mean difference. If that interval includes 0, the difference can't be told apart from chance yet: run more replicates. Also: the direct count of infections at events, and the extra infections per infection at events (above 1: the events cause more infections than happen at them; below 1: some event infections would have happened later anyway). |
| `comparison.csv` | The same values for every pair. |
| `comparison.png` | Daily new infections with and without events, and the cumulative extra infections over time (mean, with the range over pairs). |

The effect of a few events can be small compared with the normal variation between runs, so a comparison may need more replicates than a normal run. From the first event day on, the two runs of a pair use different random numbers and gradually become independent epidemics, so totals over a long run are the noisiest measure; a late event leaves less time to drift apart. When most of the population is infected either way (a high attack rate), events can barely change the total at all: they mostly change *when* people are infected. A comparison takes twice as long as a normal run with the same `n_simulations`.

## Sweeps

A sweep runs the comparison above once per value of one setting, to see how the events' effect depends on it. It's set with three things:

- `sweep_values`: the values to try, or `nothing` for no sweep (the default);
- `sweep_label`: what the values are, for the output (e.g. `"event day"`);
- `sweep`: a one-line function saying what changes for a value. It returns the settings to use; everything it doesn't return stays as set in `run.jl` and `config/user_config.jl`.

`run.jl` has two ready-made sweeps; uncomment the `sweep_values` line of the one you want:

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

`with_date(event_config, category_id, day)` puts every draw of that category on `day`, so it only works for categories with `n_draws = 1`. For a recurring category, `with_date_range(event_config, category_id, (first_day, last_day))` moves its date window instead, e.g. `sweep(start) = (; event_config = with_date_range(event_config, 2, (start, start + 20)))`. Both work for staged categories too. A sweep can change anything else `run_scenario` takes the same way, e.g. a section's size, by building the `event_config` for that value.

Every value uses the same setup seed and GEMS seed, so values differ only in the swept setting, not in luck. If `sweep` only changes `event_config` (as in sweep 1), the run without events is the same for every value, since nobody goes to any event and nothing else changes, so it runs only once: a sweep over 7 event days costs 8 runs, not 14. A sweep that changes anything else (as in sweep 2, which changes the whole epidemic) runs without events once per value. A changed transmission also moves the epidemic's peak, which changes how an event at a fixed date compares to the peak.

Results go to `Results/sweep_<n>sims_<date and time>/`:

| File | Content |
|---|---|
| `value_<value>/` | One folder per value: `with_events/` (a normal run folder), that value's `comparison.txt`, `.csv` and `.png`, and, without the shortcut, its own `without_events/`. |
| `without_events/` | The run without events shared by every value (only with the shortcut). |
| `sweep.txt` | One row per value: infections at events, extra infections (mean and 95% interval), extra among attendees and among everyone else, extra infections per event infection, and the change in peak height and day. |
| `sweep.csv` | The same values. |
| `sweep.png` | The extra infections (with the 95% interval) and the infections counted at events against the swept value, and the extra infections per event infection. |

A sweep costs (number of values + 1) runs with the shortcut, or 2 × (number of values) without, each with `n_simulations` replicates. Try it with few replicates first (see [Reproducibility and speed](reproducibility.md) for timings).

## Varying the setup per replicate

By default, all replicates of a run share one setup: the same events, the same `transmission_prob` for everyone, and the same attendees at every event. The replicates then differ only in the epidemic itself (who is infected first, who meets whom, who infects whom), so their spread shows how much the epidemic varies *given that one draw* of attendees and transmission probabilities. If that draw happens to be unusual, every replicate inherits it, and more replicates don't change that.

`vary_per_replicate` lets each replicate draw parts of the setup for itself, with the same rules:

| `vary_per_replicate` | Each replicate draws its own | Its spread then also shows |
|---|---|---|
| `[]` (default) | nothing | (only the epidemic's own randomness) |
| `[:attendance]` | who attends which event (core, loyalty, random fill, overlap, demographics as configured) | how much it matters *who* goes |
| `[:transmission]` | everyone's `transmission_prob` | how much it matters how infectious the individuals happen to be |
| `[:attendance, :transmission]` | both | both |

The events themselves (dates and sizes) are always shared, so sweeps over e.g. the event day keep the date you set. Each replicate's draw comes from its own random stream, derived from the setup seed and the replicate's number, so the same seeds give the same draws every time, and what doesn't vary is drawn exactly as without the option. GEMS's own seed per replicate doesn't change either. Comparisons and sweeps work the same way: replicate 1 with events pairs with replicate 1 without, which uses the same draws.

In the output, `attendance.jld2` holds one attendee list (and core groups) per replicate, in replicate order, and `validation.txt` checks the assignment, attendance and demographics of replicate 1's setup (every replicate's is drawn with the same rules) and uses all setups for the expected infections. The replicates are then run one by one, which costs a few extra seconds per replicate for preparing the population and assigning attendees.

## Branching

Branching runs the same simulation several times, identical up to a day you choose and different after it. Each replicate is then a **history**, and each of its runs a **future**:

```julia
branch_day = 44   # the futures differ from this day on (nothing: no branching)
n_futures = 20    # futures per history
```

Every history runs `n_futures` times. Up to day `branch_day - 1` all its futures are identical; at the end of that day, GEMS's random number generators are re-seeded, so from `branch_day` on each future goes its own way. Future 1 is not re-seeded: it is the history's own continuation, exactly as without branching. The new seeds are derived from `gems_seed`, the history and the future, so the same settings always give the same futures. A run makes `n_simulations × n_futures` GEMS runs; each future reruns the days before `branch_day`, which avoids copying a simulation in memory.

For futures of a single situation, use `n_simulations = 1`; with more, you can see whether what happens after `branch_day` depends on the history.

**Picking histories.** `histories = [7]` runs only replicate 7, exactly as it was in a run with the same `gems_seed` and settings (replicate *h* always gets the *h*-th seed GEMS derives from `gems_seed`). So to branch a replicate you saw in an earlier run, take that run's settings and `gems_seed` from its `run_info/` (see [Rerunning a run](reproducibility.md#rerunning-a-run)), set `histories` to that replicate, and choose `branch_day` and `n_futures`. Several are possible, e.g. `histories = [3, 7]`.

In the results, every future counts as a replicate; `event_results.csv` gets two extra columns, `history` (the replicate number) and `future`. Comparisons with and without events and sweeps work the same way: each history and future with events pairs with the same history and future without.
