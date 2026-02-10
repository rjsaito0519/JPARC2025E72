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
group = parser.add_mutually_exclusive_group()
group.add_argument("--bcout", action="store_true", help="Set detector to BcOut")
group.add_argument("--bcin", action="store_true", help="Set detector to BcIn")
args = parser.parse_args()

detector = "BcOut"
if args.bcin:
    detector = "BcIn"
elif args.bcout:
    detector = "BcOut"
# else: detector = "BcOut" (default)

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
            config.OUTPUT_DIR / "root" / f"run{run_num:05d}" / root_file_name,
            config.SCRATCH_DIR / f"run{run_num:05d}" / root_file_name,
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
    detectors = ["BHT", "BH2", "HTOF", "BAC", "KVC", "T1", "CVC", "SAC3", "SFV"]
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
    # C++ (BLC_TDC) outputs: runXXXXX_BLC1_TDC_Pi.root (BcIn) or runXXXXX_BLC2_TDC_Pi.root (BcOut)
    blc_det = "BLC1" if detector == "BcIn" else "BLC2" if detector == "BcOut" else detector
    
    root_file = get_root_file(args.run_num, blc_det, args.suffix, "TDC")
    if not root_file:
         root_file = get_root_file(args.run_num, detector, args.suffix, "tdc")
    if not root_file and detector != "DC":
        root_file = get_root_file(args.run_num, "DC", args.suffix, "tdc")
        
    if root_file:
        data = update_dctdc.make_dictdata(str(root_file))
        all_new_data.update(data)
        print(f"  - {detector} TDC: {color_ok()} {len(data)} entries found")

elif args.param_type == "residual":
    # Helper to apply additive update
    def prepare_additive_residuals(r_file):
        d = update_residual.make_dictdata(str(r_file))
        if not d: return {}
        
        # Get current values to ADD
        c_data = param_file.get_data_dict(key_len=cat["key_len"])
        ofs_idx = cat["start_col"] - cat["key_len"]
        
        for k in d.keys():
            if k in c_data:
                try:
                    vals = c_data[k]
                    if len(vals) > ofs_idx:
                        old_val = float(vals[ofs_idx])
                        res_val = float(d[k][0])
                        d[k] = [old_val + res_val]
                    else:
                        print(colored(f"    [Warning] Key {k}: Data row too short.", "yellow"))
                except ValueError:
                    print(colored(f"    [Warning] Key {k}: Invalid float conversion.", "yellow"))
        return d

    # C++ (BLC_residual) outputs: runXXXXX_BLC1_residual_Pi.root (BcIn) or runXXXXX_BLC2_residual_Pi.root (BcOut)
    blc_det = "BLC1" if detector == "BcIn" else "BLC2" if detector == "BcOut" else detector
    
    # 1. Search for specified detector (Result naming: BLC1_residual or BLC2_residual)
    root_file = get_root_file(args.run_num, blc_det, args.suffix, "residual")
    found_any = False
    
    if root_file:
        data = prepare_additive_residuals(root_file)
        if data:
            all_new_data.update(data)
            print(f"  - {detector} Residuals: {color_ok()} {len(data)} entries found (Additive Update)")
            found_any = True
    
    # 2. Fallbacks
    if not found_any:
        # Try "resi" suffix (Old naming or compatibility)
        root_file = get_root_file(args.run_num, detector, args.suffix, "resi")
        if root_file:
             # Double check this is not a 'decoded' file which uproot will now handle via 'tree' check
             data = prepare_additive_residuals(root_file)
             if data:
                 all_new_data.update(data)
                 print(f"  - {detector} Residuals (resi): {color_ok()} {len(data)} entries found")
                 found_any = True

    if not found_any:
        # Generic DC check
        for d_name in ["DC", "BLC1", "BLC2"]:
            for t_name in ["residual", "resi"]:
                if found_any: break
                root_file = get_root_file(args.run_num, d_name, args.suffix, t_name)
                if root_file:
                     data = prepare_additive_residuals(root_file)
                     if data:
                         all_new_data.update(data)
                         print(f"  - {d_name} {t_name}: {color_ok()} {len(data)} entries found")
                         found_any = True

# Helper to report convergence for residuals
def report_convergence(data):
    if not data: return
    
    corrs = [v[0] for v in data.values()]
    sigs = [v[1] for v in data.values() if len(v) > 1]
    
    rms = np.sqrt(np.mean(np.square(corrs)))
    max_c = np.max(np.abs(corrs))
    avg_sig = np.mean(sigs) if sigs else 0.0
    
    print(colored("\n>>> Convergence Report (Residual) <<<", "white", attrs=["bold"]))
    print(f"  - Alignment Error (RMS): {rms:.4f} mm")
    print(f"  - Maximum Correction:    {max_c:.4f} mm")
    print(f"  - Average Resolution:   {avg_sig:.4f} mm")
    
    threshold = 0.01 # 10 um
    if rms < threshold:
        print(colored(f"  - Status: [CONVERGED] (RMS < {threshold}mm)", "green", attrs=["bold"]))
    else:
        print(colored(f"  - Status: [IN PROGRESS] (Needs more iteration)", "yellow", attrs=["bold"]))
    print()

# Perform update
if all_new_data:
    # Ensure all_new_data only has one value per key for the file update 
    # (stripping extra info like sigma used for reporting)
    cleaned_data = {k: [v[0]] if isinstance(v, list) else [v] for k, v in all_new_data.items()}
    
    if args.param_type == "residual":
        report_convergence(all_new_data)
        
    success_count = param_file.update(cleaned_data, key_len=cat["key_len"], start_col=cat["start_col"])
    param_file.write()
    print(colored(f"\n[SUCCESS] Updated {success_count} entries total for Run {args.run_num} {args.suffix}", "green", attrs=["bold"]))
    print(f"File: {target_file.relative_to(config.ANALYZER_DIR)}\n")
else:
    print(colored(f"\n[SKIP] No results found for Run {args.run_num} {args.suffix} {args.param_type}.", "yellow"))
