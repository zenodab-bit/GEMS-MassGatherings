# GEMS-MassGatherings

An extension of the [GEMS](https://github.com/IMMIDD/GEMS) epidemic simulator that adds mass gathering events (festivals, concerts, sport) on top of the normal GEMS settings.

Event contacts are extra contacts: GEMS calculates each person's contacts in the normal settings (household, office, school, etc.) without considering their participation in events, and the event contacts are added on top. Events run through GEMS's global setting, which contains the entire population, so large populations are slow. Various precautions have been taken to keep it efficient, but it has not been tested on populations larger than Saalekreis.

## Setup (once per machine)

1. Install Julia 1.12 (the project was built with 1.12.5): https://julialang.org/downloads
2. Download the code and enter the folder (or download the ZIP from the GitHub page, after selecting the `multi-event-fixes` branch):
   ```
   git clone -b multi-event-fixes https://github.com/zenodab-bit/GEMS-MassGatherings.git
   cd GEMS-MassGatherings
   ```
3. Install the packages:
   ```
   julia setup.jl
   ```
   This installs the exact package versions recorded in `Manifest.toml`, GEMS included. Do not run `Pkg.update()`: it moves GEMS and the other packages to their newest versions, which may not work with this code.
4. Run a simulation, from the `GEMS-MassGatherings` folder:
   ```
   julia Concert_Project/2_Interface.jl
   ```
   Results are saved in `Results/`, inside the folder Julia was started from.

## Running a simulation

1. Configure the simulation in `Concert_Project/1_UserConfig.jl` (see [Configuring events](#configuring-events) below).
2. If needed, change the two settings in `Concert_Project/2_Interface.jl`: `n_simulations` (how many replicates to run) and `run_validation` (whether to run the checks in `7_Validation.jl`).
3. Run `2_Interface.jl`, as in step 4 of the setup.

**Reproducibility.** The same seed (set in the TOML) reproduces the same population, events and assignment every time. It does not reproduce the same epidemic outcome: GEMS's own randomness for the disease spread is separate and not tied to this seed. Re-running with the same seed gives the same setup, not the same result.

**Event dates.** Every event's date must fall within the simulation's length (`StopCriterion.limit` in the TOML). If a category's `date_range` allows a later date, `2_Interface.jl` warns but still runs. Results for that event will be inaccurate, since the simulation ends before that day is reached.

**Reading validation output.** A Z-score of "n/a" means there's only one simulation replicate, so there's no variance to compare against; it is not an error. The "Infected" row is the exception: its std comes from a binomial model, not from replicates, so it always has a Z-score. That std is a lower bound (it ignores variation in contacts and infectiousness), so treat large Z-scores there with care.

## Files

All code is in `Concert_Project/`:

| File | Purpose |
|---|---|
| `0_Helpers.jl` | Structs and helper functions used by the other files. |
| `1_UserConfig.jl` | Transmission and superspreader probabilities, and the event definitions. |
| `2_Interface.jl` | Links everything together and runs the simulations. |
| `3_Population.jl` | Creates the events and assigns the population to them. |
| `4_Contacts.jl` | Contact sampling for the whole simulation: normal settings (household, office, etc.) and events each use their own logic. |
| `5_Transmission.jl` | Transmission function for the whole simulation, with the same split as above. |
| `6_Analysis.jl` | Output: text metrics and plots. |
| `7_Validation.jl` | Checks and extra metrics to confirm the code works correctly. |
| `config_concert_covid.toml` | The standard GEMS configuration file. |

`general_rate` and `contactparameter` in the TOML must exist for GEMS's parser but are never read. The real values are set in `1_UserConfig.jl`.

## Configuring events

In `1_UserConfig.jl` you set the transmission probability for normal individuals and for superspreaders, the probability of being a superspreader, and then create categories of events. Each category has these fields:

- **`id`**: unique id used for internal processing. Must be unique, otherwise events may be overwritten.
- **`name`**: makes the outputs easier to read. Should also be unique, since some output filenames are built from it, and a collision overwrites one category's output with another's.
- **`date_range`**: the days on which events of this category can happen.
  - Events of the same category cannot happen on the same day. The code picks different dates, but the range needs enough days.
  - Different categories can land on the same date by coincidence. In that case all events of that date are filled together: core first, then loyalty, then random. So a category can never take another's core members, and its random fill can't take another's loyal attendees.
  - A person can attend events from more than one category on different days. Being core or loyal to one category doesn't exclude someone from being randomly drawn into another's events.
- **`sections`**: at least one per category. This is where the size and number of contacts are set. Multiple sections must have unique ids, and every event of the category has all its sections.
- **`n_draws`**: how many times the category recurs. The total number of events created is `n_draws` times the number of sections, not `n_draws` itself.
- **`min_age`, `max_age`**: minimum and maximum age of attendees. `age_boundaries` (see below) are shared across categories, but since the min and max age filter who is eligible first, the same boundary can cover a different real age range per category.
- **`age_weights`**: the age distribution of attendees, within the min and max ages. These are weights and don't have to sum to 1, so the real distribution at the event may differ slightly.
- **`sex_weights`**: the sex distribution of attendees. Also weights, with the same caveat.
- **`core`**: the fraction of each draw's attendees who attend every draw of the category.
  - It is calculated on the smallest draw generated (all sections together), so all core members always fit. With `section_mode = :section` it is calculated on each section's smallest size instead, so every home section fits. This gives a slightly smaller core.
  - The regular part follows the demographic weights; the superspreader part (see `min_superspreaders`) is picked uniformly at random.
  - A person can be core for at most one category.
  - With several sections, placement depends on `section_mode`.
- **`loyalty`**: the fraction of the remaining non-core spots reserved for people who already attended at least one non-core event of the category. This pool only grows and never expires. Selection follows the demographic weights; attending more times doesn't give extra weight. With several sections, placement depends on `section_mode`.
- **`min_superspreaders`**: the minimum number of superspreaders in the core, to study worst cases. More can still be picked at random later.
- **`cross_section_mean_contacts`**, **`cross_section_std_contacts`**: contacts a person has with people in other sections of the same event. The same value applies to all sections, and these contacts come on top of the within-section ones.
- **`section_mode`**: where core and loyal members go when there are several sections. `:event` (the default) places them in a random section at every draw, in proportion to section size. `:section` gives core members a fixed home section, and loyal attendees go back to their last section.

Core members and superspreaders are picked once per category, before any of its events are processed. The same core then attends every event of the category, from the first one on. Events are filled date by date: core of every event, then loyalty of every event, then random fill of every event. Anyone filled by loyalty or random fill becomes loyalty-eligible for that category's later events.

`validate_config` runs automatically once `event_config` is built. It catches common mistakes like wrong-length weight lists, inverted ranges, `core` outside [0, 1], and more `min_superspreaders` than the smallest possible draw can hold. Not everything is checked; see its comments in `0_Helpers.jl` for the full list.

Finally, `event_config` lists which categories are included in the simulation, so categories can easily be switched on or off for testing. It also holds `age_boundaries`, used to group ages the same way for all events, and `age_dist` and `sex_dist`, used for any category without its own `age_weights` or `sex_weights`.

## What happens during a run

`2_Interface.jl` loads `0_Helpers.jl` and `1_UserConfig.jl`, then the rest of the code.

**3_Population** turns the categories into actual events (each with a date and a size), builds the population, and assigns people to every event following the rules above. Nothing here needs to be configured.

Each person's transmission probability is drawn once, at the very start, from the values in `1_UserConfig.jl`, depending on whether they are a superspreader. It never changes afterward, regardless of which events that person attends. A few fields are added to every person so the rest of the code knows who attended what:

- `category_ids`, `draw_ids`, `section_ids`: which category, draw and section they attended, one entry per event attended.
- `event_dates`: the date of each event attended.
- `mean_event_contacts`, `std_event_contacts`: their within-section contact rate at each event.
- `cross_section_mean_contacts`, `cross_section_std_contacts`: the same, for contacts with people in other sections of the same event.
- `is_superspreader`, `transmission_prob`: fixed at the start, not tied to any event.

**4_Contacts** defines how many contacts a person has in every setting, both normal ones and events. For events it handles both within-section and cross-section contacts.

**5_Transmission** decides whether a contact results in an infection. Normal settings use the infecter's own `transmission_prob`. For events it also checks that the infecter and the infected attended the same event (any section).

`2_Interface.jl` then samples the events, prepares and assigns the population, and runs the GEMS simulation.

**6_Analysis** runs automatically afterwards and produces the text metrics and the plots (epidemic curves, cases by setting, per-event SEIR bars, infected boxplots), saved to the `Results` folder.

**7_Validation** runs if `run_validation` is `true`. It checks the outcomes (assignment counts, compartment totals, demographics, epidemic state) against what's expected. Not required for a normal run, but useful to catch mistakes.