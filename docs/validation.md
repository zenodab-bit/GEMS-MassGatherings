# Validation

With `run_validation = true` (in `run.jl`), every run checks its own results and writes them to `validation.txt` in the run folder (also shown on screen). It's not needed for the results themselves, but it catches mistakes in a configuration or in the code. A ✗ marks a failed check.

## What is checked

- **Assignment**: how many events people attend, that nobody attends two events on the same day (except different stages of one staged draw), and that every event got its target number of attendees.
- **Per-category attendance**: how many people attended 1, 2, … draws of each category (core and non-core), whether every core member attended every draw, and core members per section. For staged categories: each stage's attendees, and the overlap reached against the one set.
- **Results**: per event, that infected at the event never exceeds the event's size, and that the compartments add up to it.
- **Demographics**: per event, the attendees' age and sex distribution against the category's (or stage's) targets.
- **Epidemic state**: per event, the compartments and the infections at the event against expected values (below).

With `vary_per_replicate` set, the assignment, attendance and demographics checks use replicate 1's setup (every replicate's is drawn with the same rules); the epidemic-state check uses all of them.

## Reading the epidemic state

**Susceptible to Absent**: sanity checks, not formal tests. They compare the attendees with the population as a whole, assuming attendees are a random sample, which core and loyal repeat attendees aren't (they had extra exposure at earlier draws). Only very large Z-scores there deserve a closer look. Observed values are means over replicates, so each Z-score uses the spread of that mean (one run's std divided by √(number of replicates)); "n/a" means there's only one replicate, so there's no variance to compare against, not an error.

**Infected**: a real test of the event code. For each replicate, it takes who was actually present and infectious on the event day, with their own `transmission_prob`, and computes the expected infections from the rules the events are meant to follow: each infectious attendee meets a number of different people in their section or venue (drawn with its mean and std, capped at the number present), plus cross-section contacts if set, and infects each susceptible one with their own probability. It computes this directly from the distributions, without reusing the simulation's contact or transmission code, so if that code ever stops following the rules, observed and expected drift apart.

Its Z-score is the mean of the per-replicate differences (observed minus expected) over their standard error, so it needs at least 2 replicates, and more replicates give a sharper test; |Z| above about 3 is worth a closer look. On 150 replicates of the test configuration, all events stayed within ±2; with the event contacts deliberately doubled, the large events reached Z-scores of 10 and more. It doesn't test the inputs both share: the attendee lists, `transmission_prob` and the contact numbers in the configuration.

## The quick test

`test/test.jl` runs the whole pipeline on a small synthetic population (5,000 people, 60 days) with its own configuration, which uses most features: several categories and sections, both section modes, core, loyalty, cross-section contacts, and a staged category with several venues and an overlap. It also runs the comparison without events and a small sweep with per-replicate setups. Its `gems_seed` is fixed, so two test runs of the same code, with the same thread setup, give identical result files: a quick way to check that a change to the code didn't alter the results.

```
julia -t 8 test/test.jl
```

Results go to `Results/test_<n>/`. Set `test_sweep = false` in `test/test.jl` to skip the sweep (three short runs).
