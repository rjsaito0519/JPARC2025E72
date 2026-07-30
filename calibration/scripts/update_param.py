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
parser.add_argument("param_type", type=str, help="Input parameter type (hdprm, t0, hdphc, dctdc, residual, resolution)")
group = parser.add_mutually_exclusive_group()
group.add_argument("--bcout", action="store_true", help="Set detector to BcOut")
group.add_argument("--bcin", action="store_true", help="Set detector to BcIn")
group.add_argument(
    "--all",
    action="store_true",
    help="resolution only: average BLC1+BLC2 exclusive pull σ and apply same Res*=p to both",
)
args = parser.parse_args()

detector = "BcOut"
if args.bcin:
    detector = "BcIn"
elif args.bcout:
    detector = "BcOut"
# else: detector = "BcOut" (default); --all handled in resolution branch

if args.suffix not in ["K", "Pi"]:
    print("suffix should be K or Pi")
    sys.exit()

if args.all and args.param_type != "resolution":
    print(colored("[Error] --all is only supported for param_type=resolution", "red"))
    sys.exit(1)

import update_hdprm
import update_phc
import update_dctdc
import update_residual
import phc_conf
import hdprm_conf

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
    # Res column (0-based index 9): Id Name X Y Z TA RA1 RA2 L Res ...
    "resolution": {"dir": "DCGEO", "tpl": "DCGeomParam_e72_example",  "prefix": "DCGeomParam",  "key_len": 1, "start_col": 9},
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
skip_write = False
reso_meta = None

print(colored(f"\n>>> Collecting results for {args.param_type} <<<", "cyan"))

if args.param_type == "hdprm":
    detectors = ["BHT", "BH2", "HTOF", "BAC", "KVC", "T1", "CVC", "SAC3", "SFV"]
    for det in detectors:
        # HTOF HDPRM is always All (TDC-only); others use Pi/K
        suf = "All" if det == "HTOF" else args.suffix
        root_file = get_root_file(args.run_num, det, suf, "HDPRM")
        if not root_file:
            continue
            
        # Get limits
        det_lower = det.lower()
        conf_key = f"{args.run_num:05d}_{suf}_{det_lower}"
        good_range = [-np.inf, np.inf]
        if hasattr(hdprm_conf, 'limits_dict') and conf_key in hdprm_conf.limits_dict:
             good_range = hdprm_conf.limits_dict[conf_key]
             
        data = update_hdprm.make_dictdata(str(root_file), good_ch_range=good_range)
        all_new_data.update(data)
        print(f"  - {det:<6}: {color_ok()} {len(data)} entries found (Range: {good_range}, suffix={suf})")

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

elif args.param_type == "resolution":
    # exclusive Pull width p from BLC_pull -> Res' = Res * p
    import uproot

    RES_CONV_THRESHOLD = 0.05  # |p-1| below this => CONVERGED, skip Res update

    def layer_keys(blc_name):
        if blc_name == "BLC1":
            return [str(i) for i in range(1, 17)]
        if blc_name == "BLC2":
            return [str(i) for i in range(1001, 1017)]
        return None

    def find_pull_root(blc_name):
        return get_root_file(args.run_num, blc_name, args.suffix, "pull")

    def read_pull_sigma(r_file):
        """Return dict with pull_sigma, mean_chi2 (or None)."""
        try:
            file = uproot.open(str(r_file))
            if "tree" not in file:
                print(colored(f"  [Warning] 'tree' not found in {r_file}", "yellow"))
                return None
            arr = file["tree"].arrays(
                ["pull_sigma", "mean_chi2"],
                library="np",
            )
            if len(arr["pull_sigma"]) < 1:
                return None
            return {
                "pull_sigma": float(arr["pull_sigma"][0]),
                "mean_chi2": float(arr["mean_chi2"][0]),
            }
        except Exception as e:
            print(colored(f"  [Error] Reading {r_file}: {e}", "red"))
            return None

    def scale_res_for_layers(keys, p):
        """Multiply current Res by exclusive pull width p for given layer keys."""
        c_data = param_file.get_data_dict(key_len=cat["key_len"])
        res_idx = cat["start_col"] - cat["key_len"]  # index within value list
        out = {}
        sample_old = None
        sample_new = None
        for k in keys:
            if k not in c_data:
                print(colored(f"    [Warning] Key {k} not in DCGEO", "yellow"))
                continue
            vals = c_data[k]
            if len(vals) <= res_idx:
                print(colored(f"    [Warning] Key {k}: row too short for Res", "yellow"))
                continue
            try:
                old = float(vals[res_idx])
            except ValueError:
                print(colored(f"    [Warning] Key {k}: bad Res value", "yellow"))
                continue
            new = old * float(p)
            out[k] = [f"{new:.6g}"]
            if sample_old is None:
                sample_old = old
                sample_new = new
        return out, sample_old, sample_new

    def report_reso_convergence(pull_sigma, mean_chi2, sample_old, sample_new, updated):
        print(colored("\n>>> Convergence Report (Resolution / exclusive Pull) <<<", "white", attrs=["bold"]))
        print(f"  - pull σ (p):           {pull_sigma:.4f}  (target ≈ 1)")
        print(f"  - mean χ²_ν (monitor):  {mean_chi2:.4f}")
        if sample_old is not None:
            print(f"  - sample Res:           {sample_old:.4f} -> {sample_new:.4f} mm")
        dp = abs(pull_sigma - 1.0)
        if dp < RES_CONV_THRESHOLD:
            print(colored(
                f"  - Status: [CONVERGED] (|p-1|={dp:.4f} < {RES_CONV_THRESHOLD})",
                "green",
                attrs=["bold"],
            ))
            if not updated:
                print(colored("  - Res file not modified (already converged)", "green"))
        else:
            print(colored(
                f"  - Status: [IN PROGRESS] (|p-1|={dp:.4f} >= {RES_CONV_THRESHOLD})",
                "yellow",
                attrs=["bold"],
            ))
            print(colored("  - Next: re-run DST (#define FillPullExclusive 1), then reso again", "yellow"))
        print()

    found_any = False

    if args.all:
        pulls = []
        means = []
        for blc_name in ("BLC1", "BLC2"):
            root_file = find_pull_root(blc_name)
            if not root_file:
                print(colored(f"  [Error] No pull ROOT for {blc_name} (needed by --all).", "red"))
                sys.exit(1)
            info = read_pull_sigma(root_file)
            if not info:
                print(colored(f"  [Error] Cannot read pull_sigma from {root_file}", "red"))
                sys.exit(1)
            print(colored(
                f"  [INFO] {blc_name}: p={info['pull_sigma']:.4f} from {root_file.name}",
                "cyan",
            ))
            pulls.append(info["pull_sigma"])
            means.append(info["mean_chi2"])

        pull_sigma = float(np.mean(pulls)) if pulls else 0.0
        mean_chi2 = float(np.mean(means)) if means else 0.0
        print(colored(f"  [INFO] ALL mean p={pull_sigma:.4f}", "cyan"))

        reso_meta = {
            "pull_sigma": pull_sigma,
            "mean_chi2": mean_chi2,
            "sample_old": None,
            "sample_new": None,
        }

        if abs(pull_sigma - 1.0) < RES_CONV_THRESHOLD:
            skip_write = True
            found_any = True
            report_reso_convergence(pull_sigma, mean_chi2, None, None, False)
        else:
            for blc_name in ("BLC1", "BLC2"):
                keys = layer_keys(blc_name)
                data, old, new = scale_res_for_layers(keys, pull_sigma)
                all_new_data.update(data)
                if reso_meta["sample_old"] is None:
                    reso_meta["sample_old"] = old
                    reso_meta["sample_new"] = new
            print(f"  - ALL Resolution: {color_ok()} {len(all_new_data)} layers (Res *= p)")
            found_any = True
    else:
        blc_det = "BLC1" if detector == "BcIn" else "BLC2" if detector == "BcOut" else detector
        root_file = find_pull_root(blc_det)
        if not root_file:
            print(colored(f"  [Error] No pull ROOT for {blc_det}. Run BLC_pull first.", "red"))
            sys.exit(1)
        info = read_pull_sigma(root_file)
        if not info:
            print(colored(f"  [Error] Cannot read pull_sigma from {root_file}", "red"))
            sys.exit(1)

        pull_sigma = info["pull_sigma"]
        mean_chi2 = info["mean_chi2"]
        print(colored(
            f"  [INFO] {blc_det}: p={pull_sigma:.4f} mean_chi2={mean_chi2:.4f}",
            "cyan",
        ))

        reso_meta = {
            "pull_sigma": pull_sigma,
            "mean_chi2": mean_chi2,
            "sample_old": None,
            "sample_new": None,
        }

        if abs(pull_sigma - 1.0) < RES_CONV_THRESHOLD:
            skip_write = True
            found_any = True
            report_reso_convergence(pull_sigma, mean_chi2, None, None, False)
        else:
            keys = layer_keys(blc_det)
            data, old, new = scale_res_for_layers(keys, pull_sigma)
            all_new_data.update(data)
            reso_meta["sample_old"] = old
            reso_meta["sample_new"] = new
            print(f"  - {detector} Resolution: {color_ok()} {len(data)} layers (Res *= p)")
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
if args.param_type == "resolution" and skip_write:
    # Already reported CONVERGED; nothing to write
    pass
elif all_new_data:
    # Ensure all_new_data only has one value per key for the file update 
    # (stripping extra info like sigma used for reporting)
    if args.param_type == "residual":
        # report_convergence uses sigma (v[1]) if present
        report_convergence(all_new_data)
        # Stripping extra info (like sigma) for the file update
        cleaned_data = {k: [v[0]] if isinstance(v, list) else [v] for k, v in all_new_data.items()}
    else:
        cleaned_data = all_new_data
        
    success_count = param_file.update(cleaned_data, key_len=cat["key_len"], start_col=cat["start_col"])
    param_file.write()
    if args.param_type == "resolution" and reso_meta is not None:
        report_reso_convergence(
            reso_meta["pull_sigma"],
            reso_meta["mean_chi2"],
            reso_meta["sample_old"],
            reso_meta["sample_new"],
            True,
        )
    print(colored(f"\n[SUCCESS] Updated {success_count} entries total for Run {args.run_num} {args.suffix}", "green", attrs=["bold"]))
    print(f"File: {target_file.relative_to(config.ANALYZER_DIR)}\n")
else:
    print(colored(f"\n[SKIP] No results found for Run {args.run_num} {args.suffix} {args.param_type}.", "yellow"))
