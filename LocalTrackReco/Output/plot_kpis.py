import re
import matplotlib.pyplot as plt
import os


def get_event_types():
    with open(f"logs/Baseline/Baseline.log", "r") as file:
        text = file.read()
        find = re.findall(r"DownstreamTrackChecker *INFO *([0-1]\d\S*) +\:", text)
        return find


def get_thresholds(path = "./logs/Enhance"):
    filenames = [filename for filename in os.listdir(path) if os.path.isfile(os.path.join(path, filename))]
    filenames.sort()
    thresholds = [re.search(r"Enhance_(.*).log", filename).group(1) for filename in filenames]
    return thresholds


def get_value_from_file(file, regex_string):
    with open(file, "r") as file:
        text = file.read()
    value = float(re.search(regex_string, text).group(1))
    return value


def plot_pseudo_purities_for_thresholds():
    pattern = r"DownstreamTrackChecker.*\[(\d+\.\d+).*\]"
    thresholds = get_thresholds()
    pseudo_purities_for_thresholds = {}

    fraction_of_ghosts_for_baseline = get_value_from_file("logs/Baseline/Baseline.log", pattern)
    pseudo_purity_for_baseline = 1 - fraction_of_ghosts_for_baseline / 100

    for threshold in thresholds:
        fraction_of_ghosts_for_threshold = get_value_from_file(f"logs/Enhance/Enhance_{threshold}.log", pattern)
        pseudo_purities_for_thresholds[float(threshold)] = 1 - fraction_of_ghosts_for_threshold / 100

    x = pseudo_purities_for_thresholds.keys()
    y = pseudo_purities_for_thresholds.values()

    plt.plot(x, y, label="Moore Enhanced pseudo purity")
    plt.plot(x, [pseudo_purity_for_baseline for element in x], label="Moore Baseline pseudo purity", color='red', zorder=5)
    # plt.scatter(0.07, 0.792, color='red', s=20, zorder=5)
    # plt.axvline(x=0.07, color='gray', linestyle='--', linewidth=1)
    # plt.axhline(y=0.79, color='gray', linestyle='--', linewidth=1)
    # plt.text(0.2, 0.77, '(0.07, 0.79)', color='black', fontsize=10, ha='right', va='bottom')
    plt.scatter(0.34, 0.766, color='red', s=20, zorder=5)
    plt.axvline(x=0.34, color='gray', linestyle='--', linewidth=1)
    plt.axhline(y=0.766, color='gray', linestyle='--', linewidth=1)
    plt.text(0.37, 0.740, '(0.34, 0.77)', color='black', fontsize=10, ha='right', va='bottom')
    plt.title("Pseudo purities")
    plt.xlabel("threshold")
    plt.ylabel("pseudo purity")
    plt.grid(True)
    plt.legend()
    plt.savefig("pseudo_purities.png", dpi=300)
    plt.clf()


def plot_pseudo_efficiencies_for_thresholds(event_type):
    pattern = f"DownstreamTrackChecker.*{re.escape(event_type)} *: *(\d+) "
    thresholds = get_thresholds()
    pseudo_efficiencies_for_thresholds = {}

    track_number_for_baseline = get_value_from_file("logs/Baseline/Baseline.log", re.compile(pattern)) 

    for threshold in thresholds:
        track_number_for_threshold = get_value_from_file(f"logs/Enhance/Enhance_{threshold}.log", pattern)
        pseudo_efficiencies_for_thresholds[float(threshold)] = track_number_for_threshold / track_number_for_baseline

    x = pseudo_efficiencies_for_thresholds.keys()
    y = pseudo_efficiencies_for_thresholds.values()

    plt.plot(x, y, label="Moore Enhanced pseudo efficiency")
    plt.plot(x, [1 for element in x], label="Moore Baseline pseudo efficiency", color='red')
    plt.title(f"Pseudo efficiencies for {event_type}")
    plt.xlabel("threshold")
    plt.ylabel("pseudo efficiency")
    plt.grid(True)
    plt.legend()
    plt.savefig(f"pseudo_efficiencies_{event_type}.png", dpi=300)
    plt.clf()


def plot_pseudo_efficiencies_for_pseudo_purities(event_type):
    thresholds = get_thresholds()
    pattern = re.compile(f"DownstreamTrackChecker.*{re.escape(event_type)} *: *(\d+) ")
    pseudo_purities = [1 - get_value_from_file(f"logs/Enhance/Enhance_{threshold}.log", r"DownstreamTrackChecker.*\[(\d+\.\d+).*\]") / 100 for threshold in thresholds]
    track_number_for_baseline = get_value_from_file("logs/Baseline/Baseline.log", pattern) 
    pseudo_efficiencies = [get_value_from_file(f"logs/Enhance/Enhance_{threshold}.log", pattern) / track_number_for_baseline for threshold in thresholds]
    plt.plot(pseudo_efficiencies, pseudo_purities)
    plt.title(f"Pseudo efficiencies for pseudo purities for {event_type}")
    plt.xlabel("pseudo efficiency")
    plt.ylabel("pseudo purity")
    plt.grid(True)
    plt.savefig(f"pseudo_efficiencies_for_pseudo_purities_for_{event_type}.png", dpi=300)
    plt.clf()


def plot_tracks_number_for_thresholds():
    thresholds = get_thresholds()
    all_tracks_number_for_thresholds = {}
    ghost_tracks_number_for_thresholds = {}

    for threshold in thresholds:
        all_tracks_number_for_threshold = get_value_from_file(f"logs/Enhance/Enhance_{threshold}.log", "DownstreamTrackChecker.*Downstream *(\d+) ")
        all_tracks_number_for_thresholds[float(threshold)] = all_tracks_number_for_threshold

    x = all_tracks_number_for_thresholds.keys()
    y = all_tracks_number_for_thresholds.values()

    plt.plot(x, y, label="all tracks")
    plt.title("Number of tracks")
    plt.xlabel("threshold")
    plt.ylabel("number of tracks")
    plt.grid(True)

    for threshold in thresholds:
        ghost_tracks_number_for_threshold = get_value_from_file(f"logs/Enhance/Enhance_{threshold}.log", "DownstreamTrackChecker.* (\d+) ghosts")
        ghost_tracks_number_for_thresholds[float(threshold)] = ghost_tracks_number_for_threshold

    x = ghost_tracks_number_for_thresholds.keys()
    y = ghost_tracks_number_for_thresholds.values()

    plt.plot(x, y, label="ghost tracks", color="red")
    plt.legend()
    plt.savefig("number_of_tracks.png", dpi=300)
    plt.clf()


plot_pseudo_purities_for_thresholds()

# plot_pseudo_efficiencies_for_thresholds("02_UT+T_P>5GeV")
# plot_pseudo_efficiencies_for_thresholds("08_UT+T_fromBD_P>5GeV")

# plot_tracks_number_for_thresholds()

# for event_type in get_event_types():
#     plot_pseudo_efficiencies_for_thresholds(event_type)

plot_pseudo_efficiencies_for_pseudo_purities("02_UT+T_P>5GeV")
plot_pseudo_efficiencies_for_pseudo_purities("08_UT+T_fromBD_P>5GeV")