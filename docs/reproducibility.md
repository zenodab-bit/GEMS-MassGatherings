# Reproducibility and speed

## Seeds

There are two seeds:

1. **`seed` in the TOML** (`config/config_concert_covid.toml`, under `[Simulation]`) sets the setup: the events, each person's transmission probability, who attends what, and a synthetic population. The same seed gives the same setup every time. Each of these steps has its own random stream derived from `seed`, so changing one doesn't reshuffle the others: changing an event's size or `n_draws` leaves everyone's transmission probability unchanged, and changing `loyalty` or `core` leaves the events and transmission probabilities unchanged. With `vary_per_replicate`, the varying parts get one more stream per replicate, also derived from `seed`.
2. **`gems_seed` in `config/user_config.jl`** sets the epidemic itself (initial infections, contacts, infections, disease course). With `gems_seed = nothing` (the default), every run uses a new random seed, so results differ between runs. The seed used is printed and saved at the top of `multievent_analysis.txt` and in `run_info/run_info.toml`. To repeat a run exactly, set `gems_seed` to that value.

Replicates within one run always differ from each other: GEMS derives a separate seed for each one from `gems_seed` (replicate *h* always gets the *h*-th, whatever the number of replicates).

## Threads

GEMS spreads the work of each simulated day over the threads Julia was started with, so start Julia with several:

```
julia -t 8 run.jl
```

The code prints a reminder when it runs on one thread. `-t 8` is a suggestion, not a requirement: on a computer with fewer cores it still runs correctly, just no faster than with as many threads as cores (Task Manager → Performance → CPU shows the number of logical processors). Going beyond about 8 threads gains little, since only part of each day's work can be split, and a run uses those cores fully, so other programs (or other people's work on a shared server) get slower meanwhile.

The thread setup changes the random numbers GEMS draws: a `gems_seed` repeats a run exactly only with the same number of threads (and the same Julia version), whatever the number of cores. To repeat a run, use the thread count it was run with, even on a smaller computer; to just get results fast, use your own core count. Plain `julia` (without `-t`) is not the same as `julia -t 1`: since Julia 1.12 it also starts an "interactive" thread, and GEMS's results differ between the two, even though both report 1 thread. `run_info.toml` records the thread count, the interactive thread count and the Julia version of every run.

**When changing the code:** with several threads, GEMS works on different settings at the same time. A contact or transmission function must therefore not keep state that several settings share, or threads can overwrite each other's changes and silently corrupt results. The current code is safe: `NegBinContacts` keeps no state, and `EventContacts` keeps its per-day attendee lists only for the GlobalSetting, which is a single setting, so one thread at a time uses them.

## Speed

Measured on Saalekreis (434,591 people, 2 events, 120 days), 3 replicates:

| | 1 thread | 8 threads |
|---|---|---|
| Start-up: packages, population file, setup | about 25 s | about 25 s |
| First replicate (Julia compiles the code while running it) | about 58 s | about 42 s |
| Each further replicate | about 12 s | about 6.5 s |
| Saving, plots, validation (once per run) | about 20 s | about 16 s |
| Whole run | 139 s | 102 s |

So most of a short run is one-time cost, and each extra replicate is cheap. A comparison or sweep loads the population file once and reuses it for all its runs, and compiles the code only once, so it costs little more than its replicates: e.g. a sweep over 7 event days with 20 replicates (160 replicates with the shared run without events) takes about 20 minutes on 8 threads. The package itself is compiled once and cached; after a change in `src/`, the first start takes about 30 s longer.

## Rerunning a run

Every run folder has a `run_info/` folder with what's needed to repeat that run:

- copies of the files it was configured with: `config/user_config.jl`, the run script (`run.jl`, or your `my_*.jl` copy, under its own name) and the TOML (for a `test/test.jl` run: `test.jl` and the TOML);
- `run_info.toml`: the setup seed, the GEMS seed actually used, the number of replicates, the population used (file path, or the synthetic population's settings) and the settings file, the per-replicate, branching and history settings, the Julia version and thread setup, and the code version (git commit and branch);
- `uncommitted_changes.diff`, only if the code had changes not yet committed when the run started (new files never added to git are not included).

To rerun:

1. Check out the saved commit: `git checkout <git_commit>`. If there is an `uncommitted_changes.diff`, apply it: `git apply <run folder>/run_info/uncommitted_changes.diff`.
2. Copy the saved config files back: `user_config.jl` and the TOML into `config/`, the run script next to `run.jl`.
3. In `config/user_config.jl` (or `test/test.jl`), set `gems_seed` to the `gems_seed` in `run_info.toml`. The saved copy has whatever was there when the run started, usually `nothing`.
4. Make sure the population and settings files at the saved paths are the same ones, and run with the same thread setup (`julia -t <threads> ...`).

To rerun only some replicates of a run (e.g. to branch them), see `histories` in [Experiments](experiments.md#branching).
