# ===========================================================================
# config/user_config.jl
#
# The events and transmission of a simulation; how to run it is in run.jl.
# Every field is explained in src/config.jl and docs/events.md.
# ===========================================================================

# GEMS settings (pathogen, settings, simulation length, setup seed)
configfile = joinpath(@__DIR__, "config_concert_covid.toml")

# Epidemic seed (initial infections, contacts, infections, disease course):
# nothing for a new one every run (printed and saved), or a number to repeat
# a run (with the same number of threads). The setup (events, attendance,
# transmission_prob) uses the TOML's seed.
gems_seed = nothing

# Each person's transmission_prob is drawn once from this, kept within [0, 1]:
transmission_dist = Gamma(gamma_params(0.3, 0.1)...)     # by mean and std
# transmission_dist = Beta(beta_params(0.3, 0.1)...)     # by mean and std, always within [0, 1]
# transmission_dist = Dirac(0.3)                          # the same for everyone

# festival: two sections, no core or loyalty
category_1 = Category(
    id = 1,
    name = "festival",
    date_range = (40, 60),
    sections = [
        EventSection(id=1, name="seated", n_range=(300,500), mean_contacts=4.0, std_contacts=6.0),
        EventSection(id=2, name="standing", n_range=(800,1000), mean_contacts=12.0, std_contacts=18.0)
    ],
    n_draws = 0,
    min_age = 18,
    max_age = 999,
    age_weights = [0.5, 0.35, 0.15],
    sex_weights = [0.6, 0.4],
    core = 0.0,
    loyalty = 0.0,
    cross_section_mean_contacts = 1.0,
    cross_section_std_contacts = 2.0,
    section_mode = :event   # :event or :section
)

# sport: a 10% core group, one section; switched off (n_draws = 0)
category_2 = Category(
    id = 2,
    name = "sport",
    date_range = (10, 40),
    sections = [
        EventSection(id=1, name="", n_range=(300,400), mean_contacts=8.0, std_contacts=12.0)
    ],
    n_draws = 0,
    min_age = 16,
    max_age = 999,
    age_weights = [0.45, 0.40, 0.15],
    sex_weights = [0.55, 0.45],
    core = 0.1,
    loyalty = 0
)

# A staged category (docs/events.md): a match, then 10 restaurants, 30% of
# whose guests come from the match. To use it, uncomment it and add it to
# `categories` below.
 category_3 = StagedCategory(
     id = 3,
     name = "football",
     date_range = (45, 45),
     n_draws = 1,
     stages = [
         Stage(name = "match", min_age = 6,
               venues = [Venue(name = "stadium", n_range = (20000, 25000), mean_contacts = 10.0, std_contacts = 15.0)]),
         Stage(name = "restaurants", overlap = 0.3, min_age = 16,
               venues = [Venue(name = "restaurant$i", n_range = (40, 80), mean_contacts = 4.0, std_contacts = 6.0) for i in 1:2]),
     ]
 )

event_config = EventConfig(
    categories = [category_1, category_2, category_3],
    age_boundaries = [45, 65],
    age_dist = [0.211, 0.347, 0.442],
    sex_dist = [[0.420, 0.580], [0.395, 0.605], [0.426, 0.574]]
)

validate_config(event_config)

println("End user_config")