import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import numpy as np
import uproot
import os
import sys
import argparse
from pathlib import Path

# Add shared library path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config

# Apply plot settings
for k, v in config.PLOT_SETTINGS.items():
    plt.rcParams[k] = v

parser = argparse.ArgumentParser(description="Plot COBO diff")
parser.add_argument("run_num", type=int, help="Run number")
args = parser.parse_args()
run_num = args.run_num

def get_hist_data(file, key):
    if key not in file:
        print(f"Key {key} not found")
        return np.array([]), np.array([]), np.array([])
    hist_data = file[key].to_numpy()
    bin_values = hist_data[0]
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, bin_edges, bin_values

# Path format: root/runXXXXX/runXXXXX_COBO_diff.root
root_file_path = config.OUTPUT_DIR / f"root/run{run_num:05d}/run{run_num:05d}_COBO_diff.root"

if not root_file_path.exists():
    # Fallback to old path or just complain
    print(f"File not found: {root_file_path}")
    # Try flat input if it helps transition?
    flat_path = config.OUTPUT_DIR / f"root/run{run_num:05d}_COBO_diff.root"
    if flat_path.exists():
        root_file_path = flat_path
    else:
        sys.exit(1)

file = uproot.open(root_file_path)

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

center, edge, value = get_hist_data(file, f"COBO_tdc_diff_{run_num}_1")
if len(center) > 0:
    ax1.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="C0", zorder = 3)
    ax1.set_yscale("log")
    ax1.xaxis.set_major_formatter(ptick.EngFormatter())
    ax1.set_xlabel("TDC diff [ch]")

center, edge, value = get_hist_data(file, f"COBO_time_diff_{run_num}_1")
if len(center) > 0:
    ax2.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="C0", zorder = 3)
    ax2.set_yscale("log")
    ax2.set_yticklabels([])
    ax2.set_xlabel("Time diff [ns]")

plt.subplots_adjust(left = 0.1, right=0.98, top=0.92, bottom = 0.14)
img_dir = config.OUTPUT_DIR / f"img/run{run_num:05d}"
img_dir.mkdir(parents=True, exist_ok=True)
img_save_path = img_dir / f"COBO01_{run_num}.png"
plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)
# plt.show()


fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)

center, edge, value = get_hist_data(file, f"COBO_1st_hit_{run_num}")
if len(center) > 0:
    ax.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="C0", zorder = 3)
    ax.set_yscale("log")
    ax.xaxis.set_major_formatter(ptick.EngFormatter())
    ax.set_xlabel("TDC [ch]")
    ax.set_xlim(1.51E6, 1.64E6)

    value = np.array(value)
    edge  = np.array(edge)

    if np.max(value) > 0:
        threshold = 0.5 * np.max(value)
        mask = value > threshold
        if np.any(mask):
            x_min = edge[:-1][mask].min()
            x_max = edge[1:][mask].max()
            width = x_max - x_min
            print(x_min, x_max)
            print(f"Width (50% level): {width:.3e} ch")

plt.subplots_adjust(left = 0.1, right=0.98, top=0.92, bottom = 0.14)
img_save_path = img_dir / f"COBO02_{run_num}.png"
plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)
# plt.show()