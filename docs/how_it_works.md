# How it works

## The idea

The project adds mass gathering events to [GEMS](https://github.com/IMMIDD/GEMS). Event contacts come on top of the normal ones (household, office, school, ...). Events run through GEMS's GlobalSetting, which each day holds only that day's attendees, so its cost depends on event sizes, not population size. Tested on populations up to Saalekreis (about 435,000 people).

## The code

The package `MassGatherings` is in `src/`; `src/MassGatherings.jl` loads the files below in order and exports what users call. Everything else is reachable as `MassGatherings.name`.

| File | Contents |
|---|---|
| `src/utils.jl` | Distributions, age groups, plot colours, table formatting |
| `src/config.jl` | The configuration structs, `Event`, `validate_config` |
| `src/events.jl` | Event dates and sizes |
| `src/attendance.jl` | Preparing the population, assigning attendees |
| `src/contacts.jl` | Contacts; the daily GlobalSetting lists and stage passes |
| `src/transmission.jl` | Transmission |
| `src/results.jl` | Per-event counts, summaries, saving |
| `src/scenario.jl` | `run_scenario` |
| `src/comparison.jl` | `run_comparison` |
| `src/sweeps.jl` | `run_sweep`, `with_date`, `with_date_range` |
| `src/analysis.jl` | Text metrics and plots |
| `src/validation.jl` | The checks ([Validation](validation.md)) |

The package is compiled once and cached; the first start after a change in `src/` takes about 30 s longer. Changing a struct needs a new Julia session (`julia run.jl` always starts one).

## What happens during a run

`run_scenario` does these steps; comparisons, sweeps and the test all call it. Each call starts from a fresh population copy and the setup seed, so the same seeds give the same results.

**1. Events** (`events.jl`): each category becomes events with a date and size; for staged categories, one per venue.

**2. Population and attendance** (`attendance.jl`): the population passed in, or a synthetic one. Each person's `transmission_prob` is drawn once. People are assigned to events as in [Configuring events](events.md), and get these columns, which GEMS copies onto its individuals: `category_ids`, `draw_ids`, `section_ids`, `stage_ids`, `event_dates`, `mean_event_contacts`, `std_event_contacts`, `cross_section_mean_contacts`, `cross_section_std_contacts` (one entry per event attended), and `transmission_prob`.

**3. Contacts** (`contacts.jl`): an infectious attendee meets a number of different people in their section or venue (with its mean and std, capped at those present), plus cross-section contacts if set. People dead, severely ill or in hospital on the day don't attend. `daily_global_setting!` swaps each day's attendee list into the GlobalSetting (an empty one for the run without events), and stops with an error if the event columns didn't reach GEMS intact. Each further stage of a staged category gets its own pass at the end of the day; pass order doesn't matter, as people infected today only become infectious later.

**4. Transmission** (`transmission.jl`): in every setting, a contact infects with the infecter's `transmission_prob`, unless the contact has recovered.

**5. GEMS runs** the replicates: as one batch, or one by one with per-replicate setups, branching, picked histories or staged categories.

**6. Results** (`results.jl`): per event and replicate, each attendee is classified as of the event day: absent (dead, severe or hospitalised), susceptible, infectious, exposed, recovered, infected elsewhere that day, or infected at the event.

**7. Analysis** (`analysis.jl`): per event, the share of attendees present who were infected there, and the attack rate among susceptible attendees. Globally, the offspring dispersion k (lower: more superspreading), leaving out people infected too close to the end to finish infecting.

**8. Validation** (`validation.jl`), if `run_validation` is on.

## What a run folder contains

| File | Content |
|---|---|
| `multievent_analysis.txt` | Seeds, global metrics (infections, attack rate, R0, k), and per event each compartment over replicates plus the two rates; with staged categories, a per-stage summary |
| `validation.txt` | The checks (with `run_validation`) |
| `events.csv` | One row per event (per venue for staged categories, with `stage_id`), with its full `name` and its short `label` used in plots |
| `event_results.csv` | One row per replicate and event: `assigned` (people actually assigned) and every compartment; with branching or histories, `history` and `future` |
| `attendance.jld2` | `"attendees"` (event id => person ids) and `"core_groups"` (category id => person ids), one per replicate with `vary_per_replicate` |
| `*.png` | Epidemic curves (`epidemic_overview`), cases by setting, infected at each event (`infected_boxplot`) |
| `infected_by_section/` | Per section or venue, the share of susceptible attendees infected at each draw; files named `<category>_<section>` or `<category>_<stage>_<venue>` |
| `seir/` | Per event, the attendees' states on the event day; files named after the event's full name |
| `run_info/` | What's needed to rerun ([Reproducibility](reproducibility.md#rerunning-a-run)) |

Infection logs aren't saved (about 150 MB per replicate for Saalekreis).

## Populations

A population is a `.jld2` file with a table under `"data"`, a DataFrame, or `nothing` for a synthetic one. It needs `id`, `age`, `sex` (1 or 2) and the GEMS setting ids (`household`, optionally `office`, `schoolclass`, `municipality`).

The synthetic population comes from GEMS's generator (settings in `synthetic`) and the setup seed: households, offices and school classes with realistic ages, no municipalities, not census-based.

A **settings file** adds e.g. coordinates, region codes, schools and workplaces. GEMS reads it only with a population file, so the prepared population is written to a temporary file for GEMS. This makes runs noticeably slower.
