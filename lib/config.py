from pathlib import Path

# Import environment-dependent paths
try:
    from .paths import ANALYZER_DIR, PARAM_DIR, DATA_DIR, OUTPUT_DIR, SCRATCH_DIR, DECODE_DIR
except ImportError:
    raise ImportError("lib/paths.py not found. Please copy lib/paths_example.py to lib/paths.py and configure it.")

# Subdirectories
SUB_DIR = "e72"

# Detector configuration
DETECTOR_I_LIST = {
    "BHT":  1,
    "T0":   2,
    "BH2":  6,
    "BAC":  3,
    # "HTOF": 5,
    # "KVC":  6,
    # "T1":   7, 
    # "CVC":  8,
    # "SAC3": 9,
    # "SFV": 10,
}

DETECTOR_N_UD_LIST = {
    "BHT":  2,
    "T0":   2,
    "BH2":  2,
    "BAC":  1,
    "HTOF": -1,
    "KVC":  5,
    "T1":   1, 
    "CVC":  2,
    "SAC3": 1,
    "SFV":  1,
}

NUM_OF_CH = dict(
    BHT = 63,
    T0 = 5,
    BAC = 4,
    SAC3 = 1,
    KVC = 8,
    BH2 = 15,
    HTOF = 34,
    T1 = 1,
    CVC = 8,
    SFV = 1
)

# Plot settings (matplotlib)
PLOT_SETTINGS = {
    'font.family': 'Times New Roman',
    'mathtext.fontset': 'stix',
    'font.size': 22,
    'axes.linewidth': 1.0,
    'axes.grid': True,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.major.size': 10,
    'ytick.major.size': 10,
    'xtick.minor.size': 5,
    'ytick.minor.size': 5,
    'figure.subplot.left': 0.1,
    'figure.subplot.right': 0.98,
    'figure.subplot.top': 0.95,
    'figure.subplot.bottom': 0.1
}
