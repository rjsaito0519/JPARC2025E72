#!/usr/bin/env python3

# ---------------------------------------------------------------------------
import argparse
import sys
import os
import shutil
from pathlib import Path

# Add shared library path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config

# --- Color helper (No dependencies required) ---
def colored(text, color):
    codes = {"green": "32", "cyan": "36", "yellow": "33", "red": "31"}
    return f"\033[{codes.get(color, '0')}m{text}\033[0m"

parser = argparse.ArgumentParser(
    prog="make_runlist",
    usage="python3 create_runlist.py <run_num> <suffix> [--bcout | --bcin]",
    description="Tool to create analyzer .conf and .yml runlist files with smart storage management.",
    add_help=True,
)
parser.add_argument("run_num", type=int, help="Input run number")
parser.add_argument("suffix",  type=str, help="Input suffixes (e.g. Pi_hdprm K_t0)", nargs='*')

# Mutually exclusive group for detector selection
group = parser.add_mutually_exclusive_group()
group.add_argument('--bcout', action="store_true", help='Set mode to BcOut calibration')
group.add_argument('--bcin', action="store_true", help='Set mode to BcIn calibration')

args = parser.parse_args()
# ---------------------------------------------------------------------------

SUB_DIR = config.SUB_DIR

# 0. Ensure directory structures exist
(config.OUTPUT_DIR / "root" / f"run{args.run_num:05d}").mkdir(parents=True, exist_ok=True)
(config.SCRATCH_DIR / f"run{args.run_num:05d}").mkdir(parents=True, exist_ok=True)
(config.DECODE_DIR / f"run{args.run_num:05d}").mkdir(parents=True, exist_ok=True)

for param_type in ["conf", "USER", "HDPRM", "HDPHC", "DCTDC", "DCDRFT", "DCGEO"]:
    (config.PARAM_DIR / param_type / SUB_DIR).mkdir(parents=True, exist_ok=True)

# Determine Mode
# Default is Hodo
prefix = "Hodo"
mode_label = "hodo"

if args.bcout:
    prefix = "BcOut"
    mode_label = "bcout"
elif args.bcin:
    prefix = "BcIn"
    mode_label = "bcin"

set1 = set(args.suffix)
set2 = {"0", "Pi_hdprm", "K_hdprm", "Pi_t0", "K_t0", "Pi_hdphc", "K_hdphc"}
set3 = {"0", "Pi_tdc", "K_tdc", "Pi_drift", "K_drift", "Pi_resi", "K_resi"}
FINAL_SUFFIXES = ["hdphc", "resi"]

if mode_label == "hodo":
    diff = set1 - set2
    if diff or not set1:
        print(f"Error: Invalid suffixes {diff}. Use: 0, Pi/K_hdprm, Pi/K_t0, Pi/K_hdphc")
        sys.exit(1)
else:
    diff = set1 - set3
    if diff or not set1:
        print(f"Error: Invalid suffixes {diff}. Use: 0, Pi/K_tdc, Pi/K_drift, Pi/K_resi")
        sys.exit(1)

PARAM_DEFS = {
    "USER:":   {"dir": "USER",   "prefix": "UserParam_run",        "tpl": "UserParam_e72_20251104"},
    "HDPRM:":  {"dir": "HDPRM",  "prefix": "HodoParam_run",        "tpl": "HodoParam_e72_example"},
    "HDPHC:":  {"dir": "HDPHC",  "prefix": "HodoPHCParam_run",     "tpl": "HodoPHCParam_e72_example"},
    "DCTDC:":  {"dir": "DCTDC",  "prefix": "DCTdcParam_run",       "tpl": "DCTdcParam_e72_example"},
    "DCDRFT:": {"dir": "DCDRFT", "prefix": "DCDriftParam_run",     "tpl": "DCDriftParam_01928.root"},
    "DCGEO:":  {"dir": "DCGEO",  "prefix": "DCGeomParam_run",      "tpl": "DCGeomParam_e72_example"},
}

all_runs_info = []

for suffix in args.suffix:
    suffix_head = suffix.split("_")[0]
    suffix_type = suffix.split("_")[1] if "_" in suffix else suffix
    is_final = suffix_type in FINAL_SUFFIXES
    
    # 1. Update Conf File
    conf_dir = config.PARAM_DIR / "conf"
    conf_target_file = conf_dir / SUB_DIR / f"analyzer_run{args.run_num:0=5}_{mode_label}_{suffix_head}.conf"
    if not conf_target_file.exists():
        example_base = "analyzer_e72_example1.conf" if args.run_num < 2568 else "analyzer_e72_example2.conf"
        shutil.copy(conf_dir / example_base, conf_target_file)

    buf = []
    with open(conf_target_file) as f:
        for line in f:
            s_list = line.split()
            if len(s_list) > 1 and s_list[0] in PARAM_DEFS:
                p_key = s_list[0]
                p_info = PARAM_DEFS[p_key]
                filename = f"{p_info['prefix']}{args.run_num:05d}"
                if p_key != "USER:": filename += f"_{suffix_head}"
                if p_key == "DCDRFT:": filename += ".root"
                
                # Use relative path from ANALYZER_DIR for configuration consistency
                target_path_abs = config.PARAM_DIR / p_info['dir'] / SUB_DIR / filename
                if target_path_abs.exists():
                    rel_path = f"param/{p_info['dir']}/{SUB_DIR}/{filename}"
                else:
                    # Use template/default file
                    rel_path = f"param/{p_info['dir']}/{p_info['tpl']}"
                
                s_list[1] = rel_path
            buf.append(s_list)

    with open(conf_target_file, mode='w') as f:
        for l in buf:
            f.write(("\t".join(str(item) for item in l) if l else "") + "\n")

    # 2. Storage Destination - Fixed naming logic to prevent runXXXXX_DET_0_0.root
    if suffix_head == suffix_type:
        root_name = f"run{args.run_num:05d}_{prefix}_{suffix_head}.root"
    else:
        root_name = f"run{args.run_num:05d}_{prefix}_{suffix_type}_{suffix_head}.root"
        
    out_dir_abs = (config.DECODE_DIR if is_final else config.SCRATCH_DIR) / f"run{args.run_num:05d}"
    actual_root_abs = out_dir_abs / root_name
    
    # 3. Symlink at DATA_DIR root
    symlink_name = f"run{args.run_num:05d}_{prefix}.root"
    symlink_path = config.DATA_DIR / symlink_name
    if symlink_path.exists() or symlink_path.is_symlink(): symlink_path.unlink()
    try:
        os.symlink(actual_root_abs.resolve(), symlink_path)
    except Exception as e:
        print(f"Warning: Failed to create symlink: {e}")

    # Collect for runlist
    if mode_label == "hodo":
        binary = "./bin/Hodoscope"
        unit = 100000
    elif args.bcin:
        binary = "./bin/BcInTracking"
        unit = 50000
    else:
        binary = "./bin/BcOutTracking"
        unit = 50000
    all_runs_info.append({
        "label": f"run{args.run_num:05d}_{suffix}_{mode_label}:",
        "bin": binary,
        "conf": f"param/conf/{SUB_DIR}/{conf_target_file.name}",
        "data": f"rawdata/run{args.run_num:05d}.dat",
        "root": str(actual_root_abs),
        "unit": unit
    })

# 4. Update Runlist (Always recover WORKDIR and DEFAULT from myexample.yml)
runlist_dir = config.ANALYZER_DIR / "runmanager/runlist"
runlist_target_file = runlist_dir / SUB_DIR / f"{mode_label}_run{args.run_num:05d}.yml"

with open(runlist_dir / "myexample.yml") as f_tpl:
    tpl_lines = f_tpl.readlines()

with open(runlist_target_file, "w") as f_out:
    for line in tpl_lines:
        if line.strip() == "RUN:":
            break
        f_out.write(line)
    
    f_out.write("RUN:\n")
    for r in all_runs_info:
        f_out.write(f"  {r['label']}\n")
        f_out.write(f"    bin: {r['bin']}\n")
        f_out.write(f"    conf: {r['conf']}\n")
        f_out.write(f"    data: {r['data']}\n")
        f_out.write(f"    root: {r['root']}\n")
        f_out.write(f"    unit: {r['unit']}\n")

# Output Report
print(colored(f"\n[SUCCESS] Setup for Run {args.run_num} {args.suffix} ({prefix})", "green"))
print(f"  - Stable Link: {symlink_path}")

# Final verification display (Like 'cat')
print(colored("\n--- Config File Content ---", "cyan"))
sample_suffix_head = args.suffix[0].split("_")[0]
conf_file = config.PARAM_DIR / "conf" / SUB_DIR / f"analyzer_run{args.run_num:0=5}_{mode_label}_{sample_suffix_head}.conf"
print(colored(f"File: {conf_file.relative_to(config.ANALYZER_DIR)}", "yellow"))
with open(conf_file) as f: print(f.read())

print(colored("--- Runlist File Content ---", "cyan"))
print(colored(f"File: {runlist_target_file.relative_to(config.ANALYZER_DIR)}", "yellow"))
with open(runlist_target_file) as f: print(f.read())
print("-" * 60)