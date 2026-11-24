# Reproducibility and speed

## Seeds

1. **`seed` in the TOML** (`[Simulation]`) sets the setup: events, transmission probabilities, attendance and a synthetic population. Each step has its own stream from it, so changing one (e.g. an event's size) doesn't reshuffle the others.
2. **`gems_seed` in `config/user_config.jl`** sets the epidemic. `nothing` picks a new one every run; it's printed and saved in `multievent_analysis.txt` and `run_info/run_info.toml`. Set it to that value to repeat a run.

Replicate *h* always gets the *h*-th seed GEMS derives from `gems_seed`.

## Threads

```
julia -t 8 run.jl
```

More threads are faster up to about 8 (only part of each day's work can be split); with fewer cores it still works, just no faster than one thread per core. A run keeps those cores busy, so other programs, or other people's work on a shared server, slow down meanwhile.

A `gems_seed` repeats a run only with the same thread setup and Julia version, whatever the number of cores. To repeat a run, use the thread count it was run with, even on a smaller computer; to just get results fast, use your own core count. Plain `julia` isn't `julia -t 1`: since Julia 1.12 it adds an interactive thread, which changes the results. `run_info.toml` records all three.

**When changing the code**: GEMS works on several settings at once, so contact and transmission functions must not keep state shared between settings. `NegBinContacts` keeps none, and `EventContacts` only for the single GlobalSetting.

## Speed

Saalekreis (434,591 people, 2 events, 120 days), 3 replicates:

| | 1 thread | 8 threads |
|---|---|---|
| Start-up (packages, population, setup) | about 25 s | about 25 s |
| First replicate (includes compiling) | about 58 s | about 42 s |
| Each further replicate | about 12 s | about 6.5 s |
| Saving, plots, validation | about 20 s | about 16 s |
| Whole run | 139 s | 102 s |

Extra replicates are cheap. Comparisons and sweeps load the population and compile once: a sweep over 7 days with 20 replicates takes about 20 minutes on 8 threads.

## Rerunning a run

`run_info/` holds copies of the config files (`user_config.jl`, the run script, the TOML; for a test run `test.jl` and the TOML), `run_info.toml` (seeds, replicates, population, experiment settings, Julia version, threads, git commit and branch), and `uncommitted_changes.diff` if the code had uncommitted changes (new untracked files aren't included).

1. `git checkout <git_commit>`, then `git apply <run folder>/run_info/uncommitted_changes.diff` if there is one.
2. Copy the config files back.
3. Set `gems_seed` in `config/user_config.jl` (or `test/test.jl`) to the one in `run_info.toml`: the saved copy has whatever was there when the run started, usually `nothing`.
4. Use the same population files and thread setup (`julia -t <threads> ...`).

To rerun only some replicates, see `histories` in [Experiments](experiments.md#branching).
