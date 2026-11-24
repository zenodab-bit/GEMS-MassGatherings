# Validation

With `run_validation = true`, every run checks its results into `validation.txt` (also shown on screen). A ✗ marks a failed check.

## What is checked

- **Assignment**: events per person, no one at two events on a day (except stages of one staged draw), every event filled.
- **Per-category attendance**: draws attended (core and non-core), core members at every draw, core per section; for staged categories, each stage's attendees and overlap.
- **Results**: infections at an event never exceed the people assigned, and the compartments add up to them.
- **Demographics**: each event's age and sex mix against its targets.
- **Epidemic state**: compartments and infections at each event against expected values (below).

With `vary_per_replicate`, the first four use replicate 1's setup.

## Reading the epidemic state

**Susceptible to Absent**: sanity checks against the whole population, assuming attendees are a random sample, which repeat attendees aren't. Only very large Z-scores matter; "n/a" means one replicate.

**Infected**: a real test of the event code. Per replicate, it computes the expected infections from who was present and infectious, using the intended rules (contacts drawn with each section's mean and std, capped at those present, plus cross-section contacts; each infecting with the infecter's probability), without reusing the simulation code. The Z-score is the mean of the per-replicate differences over their standard error, so it needs at least 2 replicates, and few replicates give unreliable values; with many, |Z| above about 3 is worth a look. The shared inputs (attendees, `transmission_prob`, contact numbers) aren't tested.

## The quick test

`test/test.jl` runs the whole pipeline on 5,000 synthetic people for 60 days with a config using most features, then the comparison without events and a small sweep with per-replicate setups. Its `gems_seed` is fixed, so the same code and thread setup give identical result files: a quick check that a change didn't alter results.

```
julia -t 8 test/test.jl
```

Results go to `Results/test_<n>/`. `test_sweep = false` skips the sweep.
