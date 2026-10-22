# ===========================================================================
# 1_UserConfig.jl
#
# Edit this file to configure a simulation: transmission/superspreader
# parameters and event categories. This and 2_Interface.jl are the only
# files you should need to touch to run a simulation with different
# settings — everything else is internal machinery.
#
# See Section/Category/EventConfig field comments in 0_Helpers.jl for what
# each parameter below means and its constraints.
# ===========================================================================

# === Simulation ===
config = TOML.parsefile(joinpath(@__DIR__, "config_concert_covid.toml"))
seed = config["Simulation"]["seed"]
rng = Xoshiro(seed)   # same seed, same population/events/assignment every run — not the epidemic itself, GEMS seeds that separately

# Transmission — feeds gamma_params for each non-superspreader's transmission_prob
general_rate = 0.3
std_rate = 0.1

# Superspreaders — feeds gamma_params for the superspreader subgroup's transmission_prob
superspreader_prob = 0.10
superspreader_rate = 0.8
superspreader_std = 0.15

# Event definitions — two example categories below; add, remove, or edit as needed.

# festival: no core group or loyalty targeting, two sections per draw (seated/standing)
category_1 = Category(
    id = 1,
    name = "festival",
    date_range = (40, 60),
    sections = [
        Section(id=1, name="seated", n_range=(300,500), mean_contacts=4.0, std_contacts=6.0),
        Section(id=2, name="standing", n_range=(800,1000), mean_contacts=12.0, std_contacts=18.0)
    ],
    n_draws = 1,
    min_age = 18,
    max_age = 999,
    age_weights = [0.5, 0.35, 0.15],
    sex_weights = [0.6, 0.4],
    core = 0.0,
    loyalty = 0.0,
    min_superspreaders = 0,
    cross_section_mean_contacts = 1.0,
    cross_section_std_contacts = 2.0,
    section_mode = :event   # :event or :section
)

# sport: small guaranteed core group (10%) including at least 1 superspreader, recurs 5 times
category_2 = Category(
    id = 2,
    name = "sport",
    date_range = (10, 40),
    sections = [
        Section(id=1, name="", n_range=(300,400), mean_contacts=8.0, std_contacts=12.0)
    ],
    n_draws = 0,
    min_age = 16,
    max_age = 999,
    age_weights = [0.45, 0.40, 0.15],
    sex_weights = [0.55, 0.45],
    core = 0.1,
    loyalty = 0,
    min_superspreaders = 1
)

event_config = EventConfig(
    categories = [category_1, category_2],
    age_boundaries = [45, 65],
    age_dist = [0.211, 0.347, 0.442],
    sex_dist = [[0.420, 0.580], [0.395, 0.605], [0.426, 0.574]]
)

validate_config(event_config)

println("End 1_UserConfig")