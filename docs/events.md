# Configuring events

Everything here is set in `config/user_config.jl`. The structs are defined, with a comment on every field, in `src/config.jl`.

## Transmission

First, `transmission_dist`: the distribution each person's transmission probability (the chance of infecting a contact) is drawn from, once at the start. Any distribution from Distributions.jl works; values are kept within [0, 1].

```julia
transmission_dist = Gamma(gamma_params(0.3, 0.1)...)   # by mean and std (the default)
transmission_dist = Beta(beta_params(0.3, 0.1)...)     # by mean and std, always within [0, 1]
transmission_dist = Dirac(0.3)                          # the same for everyone
```

Each person's probability is used in every setting (household, office, events, ...) and never changes during a run.

## Categories

A `Category` is a kind of event that recurs (`n_draws` times), each time on its own date and with its own attendees. Example:

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

Fields:

- **`id`**: unique id used for internal processing. `validate_config` stops with an error if two categories share an id.
- **`name`**: makes the outputs easier to read. Should also be unique, since some output filenames are built from it, and a collision overwrites one category's output with another's.
- **`date_range`**: the days on which events of this category can happen.
  - Events of the same category cannot happen on the same day. The code picks different dates, but the range needs enough days.
  - Different categories can land on the same date by coincidence. In that case all events of that date are filled together: core first, then loyalty, then random. So a category can never take another's core members, and its random fill can't take another's loyal attendees.
  - A person can attend events from more than one category on different days. Being core or loyal to one category doesn't exclude someone from being randomly drawn into another's events.
- **`sections`**: at least one `EventSection` per category, with its own size (`n_range`, drawn anew at every draw) and contacts within it (`mean_contacts`, `std_contacts`; a negative binomial, or Poisson when std² ≤ mean). Multiple sections must have unique ids, and every event of the category has all its sections.
- **`n_draws`**: how many times the category recurs. The total number of events created is `n_draws` times the number of sections, not `n_draws` itself. A negative value is an error; 0 gives a warning and creates no events for that category.
- **`min_age`, `max_age`**: minimum and maximum age of attendees. `age_boundaries` (below) are shared across categories, but since the min and max age filter who is eligible first, the same boundary can cover a different real age range per category.
- **`age_weights`**: the age distribution of attendees, within the min and max ages. These are weights and don't have to sum to 1, so the real distribution at the event may differ slightly.
- **`sex_weights`**: the sex distribution of attendees. Also weights, with the same caveat.
- **`core`**: the fraction of each draw's attendees who attend every draw of the category.
  - It is calculated on the smallest draw generated (all sections together), so all core members always fit. With `section_mode = :section` it is calculated on each section's smallest size instead, so every home section fits. This gives a slightly smaller core.
  - Core members are picked once, before any draw, following the demographic weights. A person can be core for at most one category.
- **`loyalty`**: the fraction of the spots left after core reserved for people who already attended at least one earlier draw of the category as non-core. This pool only grows and never expires. Selection follows the demographic weights; attending more times doesn't give extra weight. At the first draw there's no pool yet, so loyalty has no effect there.
- **`cross_section_mean_contacts`**, **`cross_section_std_contacts`**: contacts a person has with people in other sections of the same draw, on top of the within-section ones. The same value applies to all sections; 0 keeps sections apart.
- **`section_mode`**: where core and loyal members go when there are several sections. `:event` (the default) places them in a random section at every draw, in proportion to section size. `:section` gives core members a fixed home section, and loyal attendees go back to their last section.

Events are filled date by date: core of every event, then loyalty of every event, then random fill of every event. Anyone filled by loyalty or random fill becomes loyalty-eligible for that category's later draws. If not enough eligible people are left on a day to fill an event (e.g. small population, strict age limits, several big events on the same day), the event is filled as far as possible and a warning is printed.

Example of core and loyalty: draws of 1,000 people with `core = 0.3` and `loyalty = 0.3` have 300 core members at every draw; from the second draw on, 30% of the 700 remaining spots (210) go to people who were there before, and 490 are random. So about 51% of each later draw has been there before.

## Staged categories

A `StagedCategory` is a gathering made of several **stages on the same day**, where people can attend more than one stage: e.g. restaurants that feed a football match. Each stage has one or more **venues** (several restaurants, or one or more stadiums); a venue is a single group, where the people present meet each other, with no contacts to other venues and no internal sections.

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

(This example is in `config/user_config.jl`, commented out.)

- **`StagedCategory`**: `id`, `name`, `date_range` and `n_draws` as for a `Category`, and its `stages` in order. All stages of a draw are on the same day. There is no core or loyalty: every draw draws its attendees afresh.
- **`Stage`**: `name`, `venues`, `overlap`, and its own `min_age`, `max_age`, `age_weights` and `sex_weights` (as for a `Category`; without weights, `age_dist` and `sex_dist` are used).
- **`overlap`**: the share of the stage's attendees taken from the **previous** stage's attendees of the same draw; it must be 0 for the first stage. The order of the stages decides the direction: in the example, 30% of the restaurant guests come from the match and 70% didn't go to it. With the restaurants first and the match second, the match's overlap would say how many match-goers come from the restaurants instead. If the previous stage hasn't enough eligible people, the overlap is filled as far as possible and a warning is printed.
- **`Venue`**: `name`, `n_range` (its size, drawn anew at every draw), and the contacts within it (`mean_contacts`, `std_contacts`, as for a section).

How the attendees are drawn, per draw and stage: first the overlap from the previous stage, then the rest from people at no earlier stage of the draw and at no other event that day, following the stage's demographic weights; then everyone is spread over the stage's venues in proportion to their sizes, anew at every draw. Nobody is at two venues of the same stage. On a day shared with normal categories, staged categories are filled after them.

During the simulation each stage gets its own pass through GEMS's GlobalSetting (see [How it works](how_it_works.md)), so the results can tell the stages apart: the per-event results are per venue, `multievent_analysis.txt` gets a summary of the infections per stage, and the validation checks the overlap reached and the expected infections per venue. Two people who meet both at a restaurant and at the stadium have two chances of infection, one per stage.

## The event configuration

Finally, `event_config` lists which categories are included in the simulation, so categories can easily be switched on or off:

```julia
event_config = EventConfig(
    categories = [category_1, category_2],
    age_boundaries = [45, 65],
    age_dist = [0.211, 0.347, 0.442],
    sex_dist = [[0.420, 0.580], [0.395, 0.605], [0.426, 0.574]]
)
validate_config(event_config)
```

- `age_boundaries` group ages the same way for all events (here: up to 45, 46–65, over 65).
- `age_dist` and `sex_dist` are used for any category or stage without its own `age_weights` or `sex_weights`.
- `validate_config` catches common mistakes: wrong-length weight lists, inverted ranges, `core` or `overlap` outside [0, 1], an overlap on a first stage, negative `n_draws`, duplicate category or section ids. Not everything is checked; see its comments in `src/config.jl`.

## Event dates

GEMS simulates days 0 to `StopCriterion.limit` − 1 (in the TOML), so with `limit = 120` the last simulated day is 119. Every event's date must be on or before that day. If a category's `date_range` allows a later date, the run warns but still goes ahead; results for that event will be inaccurate, since the simulation ends before that day is reached.

## The GEMS configuration (TOML)

`config/config_concert_covid.toml` is a standard GEMS configuration file: pathogen, settings, simulation length, start condition and the setup seed. In it, the transmission function's `parameters` table is empty on purpose: GEMS requires the table, but transmission uses each person's own probability, set in `config/user_config.jl`. Likewise, the GlobalSetting's contact method has no parameters, since event contact rates are set per section and venue in `config/user_config.jl`.
