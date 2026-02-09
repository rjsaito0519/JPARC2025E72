#!/usr/bin/env python3

import argparse
import sys
import os
import shutil
import numpy as np
from termcolor import colored
from pathlib import Path

# Add shared library path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config
from lib.param import ParamFile

parser = argparse.ArgumentParser(
    prog="update_param",
    description="Update parameter files from ROOT analysis results",
    add_help=True,
)
parser.add_argument("run_num", type=int, help="Input run number")
parser.add_argument("suffix", type=str, help="Input suffix (K or Pi)")
parser.add_argument("param_type", type=str, help="Input parameter type (hdprm, t0, hdphc, dctdc, residual)")
args = parser.parse_args()

if args.suffix not in ["K", "Pi"]:
    print("suffix should be K or Pi")
    sys.exit()

import update_hdprm
import update_phc
import update_dctdc
import update_residual
import phc_conf

def get_root_file(run_num, det, suffix, type_str):
    """Search for ROOT file in scratch, decode, or local output dirs."""
    patterns = []
    
    # Construct base file name components
    parts = [f"run{run_num:0=5}", det]
    if type_str: parts.append(type_str)
    if suffix: parts.append(suffix)
    
    # 1. Preferred pattern (components joined by underscore)
    patterns.append("_".join(parts) + ".root")
    
    # 2. Alternative patterns for robustness
    if type_str and suffix:
        # Try swapping type and suffix if needed
        alt_parts = [f"run{run_num:0=5}", det, suffix, type_str]
        patterns.append("_".join(alt_parts) + ".root")
    
    # 3. Simplest pattern
    if suffix:
        patterns.append(f"run{run_num:0=5}_{det}_{suffix}.root")
    else:
        patterns.append(f"run{run_num:0=5}_{det}.root")

    for root_file_name in patterns:
        search_paths = [
            config.SCRATCH_DIR / f"run{run_num:05d}" / root_file_name,
            config.DECODE_DIR / f"run{run_num:05d}" / root_file_name,
            config.OUTPUT_DIR / "root" / f"run{run_num:05d}" / root_file_name
        ]
        for root_file in search_paths:
            if root_file.exists():
                return root_file
    return None

def color_ok():
    # Force green [OK] for visibility even if redirected
    return f"[{colored('OK', 'green', attrs=['bold'])}]"

# Dispatch configurations and Key Lengths
CATEGORIES = {
    "hdprm":    {"dir": "HDPRM", "tpl": "HodoParam_e72_example",      "prefix": "HodoParam",    "key_len": 5, "start_col": 5},
    "t0":       {"dir": "HDPRM", "tpl": "HodoParam_e72_example",      "prefix": "HodoParam",    "key_len": 5, "start_col": 5},
    "hdphc":    {"dir": "HDPHC", "tpl": "HodoPHCParam_e72_example",   "prefix": "HodoPHCParam", "key_len": 4, "start_col": 4},
    "dctdc":    {"dir": "DCTDC", "tpl": "DCTdcParam_e72_example",     "prefix": "DCTdcParam",   "key_len": 3, "start_col": 3},
    "residual": {"dir": "DCGEO", "tpl": "DCGeomParam_e72_example",    "prefix": "DCGeomParam",  "key_len": 1, "start_col": 12},
}

if args.param_type not in CATEGORIES:
    print(colored(f"Error: Unknown param_type {args.param_type}", "red"))
    sys.exit(1)

cat = CATEGORIES[args.param_type]
SUB_DIR = config.SUB_DIR
target_dir = config.PARAM_DIR / cat["dir"] / SUB_DIR
target_dir.mkdir(parents=True, exist_ok=True)
target_file = target_dir / f"{cat['prefix']}_run{args.run_num:0=5}_{args.suffix}"
template_file = config.PARAM_DIR / cat["dir"] / cat["tpl"]

param_file = ParamFile.create_from_template(template_file, target_file)

# Collect new data
all_new_data = {}

print(colored(f"\n>>> Collecting results for {args.param_type} <<<", "cyan"))

if args.param_type == "hdprm":
    detectors = ["BHT", "BH2", "BAC", "KVC", "T1", "CVC", "SAC3", "SFV"]
    for det in detectors:
        root_file = get_root_file(args.run_num, det, args.suffix, "HDPRM")
        if not root_file:
            continue
            
        # Get limits
        det_lower = det.lower()
        conf_key = f"{args.run_num:05d}_{args.suffix}_{det_lower}"
        good_range = [-np.inf, np.inf]
        if hasattr(phc_conf, 'limits_dict') and conf_key in phc_conf.limits_dict:
             good_range = phc_conf.limits_dict[conf_key]
             
        data = update_hdprm.make_dictdata(str(root_file), good_ch_range=good_range)
        all_new_data.update(data)
        print(f"  - {det:<6}: {color_ok()} {len(data)} entries found (Range: {good_range})")

elif args.param_type == "t0":
    root_file = get_root_file(args.run_num, "T0_Offset", args.suffix, "")
    if root_file:
        data = update_hdprm.make_dictdata(str(root_file), is_t0_offset=True)
        all_new_data.update(data)
        print(f"  - T0_Offset : {color_ok()} {len(data)} entries found")

elif args.param_type == "hdphc":
    detectors = ["BHT", "BH2", "HTOF", "T1", "CVC"]
    for det in detectors:
        root_file = get_root_file(args.run_num, det, args.suffix, "PHC")
        if not root_file:
            continue
            
        # Get limits
        det_lower = det.lower()
        conf_key = f"{args.run_num:05d}_{args.suffix}_{det_lower}"
        good_range = [-np.inf, np.inf]
        if hasattr(phc_conf, 'limits_dict') and conf_key in phc_conf.limits_dict:
             good_range = phc_conf.limits_dict[conf_key]

        data = update_phc.make_dictdata(str(root_file), good_ch_range=good_range)
        all_new_data.update(data)
        print(f"  - {det:<6}: {color_ok()} {len(data)} entries found (Range: {good_range})")

elif args.param_type == "dctdc":
    root_file = get_root_file(args.run_num, "DC", args.suffix, "tdc")
    if not root_file:
        root_file = get_root_file(args.run_num, "BLC2", args.suffix, "TDC")
    if root_file:
        data = update_dctdc.make_dictdata(str(root_file))
        all_new_data.update(data)
        print(f"  - DC/BLC TDC: {color_ok()} {len(data)} entries found")

elif args.param_type == "residual":
    root_file = get_root_file(args.run_num, "DC", args.suffix, "resi")
    if root_file:
        data = update_residual.make_dictdata(str(root_file))
        all_new_data.update(data)
        print(f"  - DC Residuals: {color_ok()} {len(data)} entries found")
    else:
        for det in ["BLC1", "BLC2"]:
            root_file = get_root_file(args.run_num, det, args.suffix, "residual")
            if root_file:
                data = update_residual.make_dictdata(str(root_file))
                all_new_data.update(data)
                print(f"  - {det:<6}: {color_ok()} {len(data)} entries found")

# Perform update
if all_new_data:
    success_count = param_file.update(all_new_data, key_len=cat["key_len"], start_col=cat["start_col"])
    param_file.write()
    print(colored(f"\n[SUCCESS] Updated {success_count} entries total for Run {args.run_num} {args.suffix}", "green", attrs=["bold"]))
    print(f"File: {target_file.relative_to(config.ANALYZER_DIR)}\n")
else:
    print(colored(f"\n[SKIP] No results found for Run {args.run_num} {args.suffix} {args.param_type}.", "yellow"))
