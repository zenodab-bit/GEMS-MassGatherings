# Configuring events

Everything here is set in `config/user_config.jl`; the structs, with a comment on every field, are in `src/config.jl`.

## Transmission

`transmission_dist` is the distribution each person's transmission probability is drawn from, once. Any Distributions.jl distribution works; values are kept within [0, 1]. It applies in every setting and never changes during a run.

```julia
transmission_dist = Gamma(gamma_params(0.3, 0.1)...)   # by mean and std (the default)
transmission_dist = Beta(beta_params(0.3, 0.1)...)     # by mean and std, always within [0, 1]
transmission_dist = Dirac(0.3)                          # the same for everyone
```

## Categories

A `Category` is a kind of event that recurs `n_draws` times, each on its own date with its own attendees:

```julia
category_1 = Category(
    id = 1,
    name = "festival",
    date_range = (40, 60),
    sections = [
        EventSection(id=1, name="seated",   n_range=(300,500),  mean_contacts=4.0,  std_contacts=6.0),
        EventSection(id=2, name="standing", n_range=(800,1000), mean_contacts=12.0, std_contacts=18.0)
    ],
    n_draws = 1,
    min_age = 18,
    age_weights = [0.5, 0.35, 0.15],
    sex_weights = [0.6, 0.4],
    cross_section_mean_contacts = 1.0,
    cross_section_std_contacts = 2.0,
)
```

- **`id`**: unique.
- **`name`**: event names are built from it and the draw and section names, and plot files are named after them, so no two events may share a name. In plots, events are shown by a short label: the section's (or for staged categories the venue's) own name, with the draw number when there are several draws, e.g. `seated_2`; the category name if the section has none; the full name if two events would share a label.
- **`date_range`**: the days its draws can happen on, one draw per day, so it needs at least `n_draws` days. Categories can share a date; their events are then filled together (core, then loyalty, then random fill), so no category takes another's core or loyal attendees. People can attend different categories on different days.
- **`sections`**: one or more, each with its own size (`n_range`, drawn at every draw) and contacts (`mean_contacts`, `std_contacts`: negative binomial, or Poisson when std² ≤ mean). Every draw has all sections.
- **`n_draws`**: how often it recurs; 0 switches it off (with a warning).
- **`min_age`, `max_age`**: who is eligible.
- **`age_weights`, `sex_weights`**: the attendees' demographics, as weights (they needn't sum to 1). A weight of 0 leaves a group out.
- **`core`**: the share of each draw who attend every draw, computed on the smallest draw (with `section_mode = :section`, on each section's smallest size) so it always fits. Picked once; nobody is core for two categories.
- **`loyalty`**: the share of the spots left after core kept for people who attended an earlier draw as non-core. The pool only grows. No effect at the first draw.
- **`cross_section_mean_contacts`, `cross_section_std_contacts`**: contacts with other sections of the same draw; 0 keeps sections apart.
- **`section_mode`**: with several sections, `:event` (default) puts core and loyal members in a random section each draw; `:section` gives core a fixed section and sends loyal attendees back to their last one.

Example: draws of 1,000 with `core = 0.3` and `loyalty = 0.3` have 300 core members; from the second draw, 210 of the 700 other spots go to returning attendees and 490 are random.

If too few eligible people are left on a day, an event is filled as far as possible with a warning.

## Staged categories

A `StagedCategory` is several **stages on the same day**, where people can attend more than one stage, e.g. restaurants and a football match. Each stage has **venues**: single groups with no contacts between them.

```julia
category_3 = StagedCategory(
    id = 3,
    name = "football",
    date_range = (20, 50),
    n_draws = 2,
    stages = [
        Stage(name = "match", min_age = 6,
              venues = [Venue(name = "stadium", n_range = (20000, 25000), mean_contacts = 10.0, std_contacts = 15.0)]),
        Stage(name = "restaurants", overlap = 0.3, min_age = 16,
              venues = [Venue(name = "restaurant$i", n_range = (40, 80), mean_contacts = 4.0, std_contacts = 6.0) for i in 1:10]),
    ]
)
```

- **`StagedCategory`**: `id`, `name`, `date_range` and `n_draws` as for a `Category`, and `stages` in order. No core or loyalty.
- **`Stage`**: `name`, `venues`, `overlap`, and its own age limits and weights.
- **`overlap`**: the share of the stage's attendees taken from the previous stage's (0 for the first). The stage order sets the direction: here 30% of the restaurant guests were at the match and 70% weren't. With the restaurants first and the match second, the match's overlap would instead say how many match-goers came from the restaurants. If the previous stage has too few eligible people, the overlap is filled as far as possible with a warning.
- **`Venue`**: `name`, `n_range` and contacts, as for a section.

Per stage, the overlap comes first, then people at no earlier stage and no other event that day; then they are spread over the venues by size. Staged categories are filled after normal ones on the same day. Each stage gets its own pass through GEMS's GlobalSetting, so results are per venue and stage. Two people who meet both at a restaurant and at the stadium have two chances of infection, one per stage.

## The event configuration

```julia
event_config = EventConfig(
    categories = [category_1, category_2],
    age_boundaries = [45, 65],
    age_dist = [0.211, 0.347, 0.442],
    sex_dist = [[0.420, 0.580], [0.395, 0.605], [0.426, 0.574]]
)
validate_config(event_config)
```

- `categories`: the ones included, so they're easy to switch on and off.
- `age_boundaries`: age groups for all events (here up to 45, 46–65, over 65), increasing.
- `age_dist`, `sex_dist`: for categories or stages without their own weights.
- `validate_config` stops on common mistakes (wrong lengths, inverted or negative ranges, shares outside [0, 1], duplicate ids or event names); `run_scenario` runs it again before every run.

## Event dates

GEMS simulates days 0 to `StopCriterion.limit` − 1. Events after that warn and give inaccurate results.

## The GEMS configuration (TOML)

`config/config_concert_covid.toml` is a standard GEMS configuration: pathogen, settings, simulation length, start condition and setup seed. The transmission function's `parameters` and the GlobalSetting's contact method are empty on purpose: both are set in `config/user_config.jl`.
