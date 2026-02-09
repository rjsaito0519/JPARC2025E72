#!/usr/bin/env python3

# ---------------------------------------------------------------------------
import argparse
import sys
import os
import shutil
from pathlib import Path

# Immediate feedback
print("\033[36mStarting up set_tdc_window tool...\033[0m")

import numpy as np
import uproot
import matplotlib.pyplot as plt
from termcolor import colored

# Add shared library path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config
from lib.param import ParamFile

parser = argparse.ArgumentParser(
    prog="set_tdc_window",
    usage="python3 set_tdc_window.py <rootfile_path> <counter_name>",
    description="Interactive TDC window setting tool",
    add_help=True,
)
parser.add_argument("rootfile_path", type=str, help="Input rootfile path")
parser.add_argument("counter_name", type=str, help="Input counter name")
args = parser.parse_args()
# ---------------------------------------------------------------------------

filename = os.path.basename(args.rootfile_path)
index = filename.find("run")
if index == -1:
    print("Error: Could not find 'run' in filename to parse run number.")
    run_num = 0
else:
    run_num = int(filename[index+3:index+8])

SUB_DIR = config.SUB_DIR
user_param_dir = config.PARAM_DIR / "USER" / SUB_DIR
user_param_dir.mkdir(parents=True, exist_ok=True)
target_file = user_param_dir / f"UserParam_run{run_num:05d}"

if not target_file.exists():
    shutil.copy(config.PARAM_DIR / "USER/UserParam_e72_example", target_file)

# Apply plot settings from config
for k, v in config.PLOT_SETTINGS.items():
    plt.rcParams[k] = v

# matplotlibのグラフを動的にする
# ---------------------------------------------------------------------------
fill_patch = None
vline_patch = None
bht_mode = False
selected_points = []

def on_click(event):
    global fill_patch
    global vline_patch
    global bht_mode

    if event.button == 3: # Right click
        x = event.xdata
        if x is None: return
        print(f'Right click at x={x:.2f}')
        selected_points.append(x)
        
        if len(selected_points) > 1:
            if fill_patch: fill_patch.remove()
            if vline_patch: vline_patch.remove()
            
            if bht_mode and len(selected_points) >= 3:
                x1, x2, x3 = selected_points[-3], selected_points[-2], selected_points[-1]
                x_min, x_max = sorted([x1, x2])
                y_min, y_max = ax.get_ylim()
                fill_patch = ax.fill_betweenx([y_min, y_max], x_min, x_max, color='gray', alpha=0.5)
                vline_patch = ax.vlines(x3, y_min, y_max, color='red', ls="dashed")
            elif not bht_mode:
                x1, x2 = selected_points[-2], selected_points[-1]
                x_min, x_max = sorted([x1, x2])
                y_min, y_max = ax.get_ylim()
                fill_patch = ax.fill_betweenx([y_min, y_max], x_min, x_max, color='gray', alpha=0.5)
            
            plt.draw()

def on_key(event):
    global bht_mode          
    if event.key == "w":
        if (bht_mode and len(selected_points) >= 3) or (not bht_mode and len(selected_points) >= 2):
            tdc_label = f"{args.counter_name}_TDC"
            updated_params = {}
            
            if bht_mode:
                x1, x2, x3 = selected_points[-3], selected_points[-2], selected_points[-1]
                xmin, xmax = sorted([x1, x2])
                p_sep = int(x3 + 0.5)
                p_min, p_max = int(xmin + 0.5), int(xmax + 0.5)
                updated_params[tdc_label] = [p_min, p_max]
                updated_params[f"{tdc_label}_Pi"] = [p_min, p_sep]
                updated_params[f"{tdc_label}_K"] = [p_sep, p_max]
            else:
                x1, x2 = selected_points[-2], selected_points[-1]
                p_min, p_max = sorted([int(x1 + 0.5), int(x2 + 0.5)])
                updated_params[tdc_label] = [p_min, p_max]

            # Use ParamFile for robust update
            pfile = ParamFile(target_file)
            pfile.update(updated_params, key_len=1)
            pfile.write()
            
            print(f"==========\nUpdated {target_file.name}:")
            for k, v in updated_params.items():
                print(f"  {k:20}: {v}")
            print("==========")
                
    elif event.key == "1":
        print("bht mode ON")
        bht_mode = True
    elif event.key == "2":
        print("bht mode OFF")
        bht_mode = False
# ---------------------------------------------------------------------------
print(f"Loading histograms for {args.counter_name}...", end="", flush=True)
file = uproot.open(args.rootfile_path)

# Handle different detector types
if args.counter_name in ["BHT", "T0", "T1", "CVC", "BH2", "HTOF", "SAC3", "SFV"]:
    hist_key_prefix = f"{args.counter_name}_TDC_seg"
    num_ch = config.NUM_OF_CH[args.counter_name]
    bin_values_u = None
    bin_values_d = None
    bin_edges = None

    for i in range(num_ch):
        hu_key = f"{hist_key_prefix}{i}U"
        hd_key = f"{hist_key_prefix}{i}D"
        if hu_key in file:
            hu = file[hu_key].to_numpy()
            if bin_values_u is None:
                bin_values_u = np.copy(hu[0])
                bin_edges = hu[1]
            else:
                bin_values_u += hu[0]
        if hd_key in file:
            hd = file[hd_key].to_numpy()
            if bin_values_d is None:
                bin_values_d = np.copy(hd[0])
            else:
                bin_values_d += hd[0]
    
    print(" Done.")
    if bin_values_u is None:
        print(colored(f"\nError: No histograms found for {args.counter_name} in {args.rootfile_path}", "red"))
        sys.exit(1)
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)

    ax.hist(bin_centers, bins=bin_edges, weights=bin_values_u, histtype='step', color="C1", label="UP")
    ax.hist(bin_centers, bins=bin_edges, weights=bin_values_d, histtype='step', color="C0", label="DOWN")
    ax.set_title(f"TDC Window Setup: {args.counter_name} (Run {run_num})")
    ax.set_xlabel("TDC [channel]")
    ax.set_yscale("log")

elif args.counter_name in ["BAC", "SAC3"]:
    num_ch = config.NUM_OF_CH[args.counter_name]
    hist_sum_key = f"{args.counter_name}_TDC_seg{num_ch}U"
    
    bin_values = None
    bin_edges = None

    if hist_sum_key in file:
        h = file[hist_sum_key].to_numpy()
        bin_values = h[0]
        bin_edges = h[1]
    else:
        for i in range(num_ch):
            h_key = f"{args.counter_name}_TDC_seg{i}U"
            if h_key in file:
                h = file[h_key].to_numpy()
                if bin_values is None:
                    bin_values = np.copy(h[0])
                    bin_edges = h[1]
                else:
                    bin_values += h[0]
    
    print(" Done.")
    if bin_values is None:
        print(colored(f"Error: No histograms found for {args.counter_name} in {args.rootfile_path}", "red"))
        sys.exit(1)

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.hist(bin_centers, bins=bin_edges, weights=bin_values, histtype='step')
    ax.set_title(f"TDC Window Setup: {args.counter_name} (Run {run_num})")
    ax.set_yscale("log")
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)

elif args.counter_name == "KVC":
    num_ch = config.NUM_OF_CH["KVC"]
    bin_values = None
    bin_edges = None
    for i in range(num_ch):
        h_key = f"KVC_TDC_seg{i}S"
        if h_key in file:
            h = file[h_key].to_numpy()
            if bin_values is None:
                bin_values = np.copy(h[0])
                bin_edges = h[1]
            else:
                bin_values += h[0]
            
    print(" Done.")
    if bin_values is None:
        print(colored(f"Error: No KVC histograms found in {args.rootfile_path}", "red"))
        sys.exit(1)

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.hist(bin_centers, bins=bin_edges, weights=bin_values, histtype='step')
    ax.set_title(f"TDC Window Setup: KVC (Run {run_num})")
    ax.set_yscale("log")
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)

# Show current guidelines
if target_file.exists():
    pfile = ParamFile(target_file)
    cur_data = pfile.get_data_dict()
    tdc_label = f"{args.counter_name}_TDC"
    y_min, y_max = ax.get_ylim()

    if tdc_label in cur_data:
        vals = cur_data[tdc_label]
        x_min, x_max = float(vals[0]), float(vals[1])
        ax.fill_betweenx([y_min, y_max], x_min, x_max, color='gray', alpha=0.15, label="Current Window")
        
        if args.counter_name == "BHT":
            pi_label = f"{tdc_label}_Pi"
            if pi_label in cur_data:
                x_sep = float(cur_data[pi_label][1])
                ax.vlines(x_sep, y_min, y_max, color='red', ls="dotted", alpha=0.5, label="Current Pi/K Sep")
    ax.legend()

plt.show()
