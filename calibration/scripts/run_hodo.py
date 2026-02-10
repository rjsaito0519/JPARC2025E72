#!/usr/bin/env python3

import argparse
import sys
import os
import subprocess
from pathlib import Path
from termcolor import colored

# Add shared library path to find config
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config

def run_command(cmd, exit_on_error=True):
    print(colored(f"Running: {cmd}", "cyan"))
    ret = subprocess.call(cmd, shell=True)
    if ret != 0:
        print(colored(f"[Error] Command failed with exit code {ret}", "red"))
        if exit_on_error:
            sys.exit(ret)
    return ret

def main():
    parser = argparse.ArgumentParser(
        description="Automate Hodoscope Calibration for HDPRM, T0 Offset, and HDPHC."
    )
    parser.add_argument("run_num", type=int, help="Run Number")
    parser.add_argument("mode", type=str, choices=["hdprm", "t0", "hdphc"], help="Calibration Mode")
    parser.add_argument('--kaon', action="store_true", help='Use Kaon (K) suffix instead of Pion (Pi)')
    
    args = parser.parse_args()

    run_num = args.run_num
    mode = args.mode
    suffix = "K" if args.kaon else "Pi"

    # 1. Locate the Input Root File (Symbolic link in DATA_DIR)
    # create_runlist.py creates runXXXXX_Hodo.root
    
    root_filename = f"run{run_num:05d}_Hodo.root"
    input_root_file = config.DATA_DIR / root_filename
    
    if not input_root_file.exists():
         print(colored(f"[Error] Symlink/File not found: {input_root_file}", "red"))
         print(f"Ensure create_runlist.py was run for Hodo.")
         sys.exit(1)
         
    print(colored(f"[INFO] Using Input File: {input_root_file}", "green"))
    print(colored(f"[INFO] Mode: {mode} (Suffix: {suffix})", "green"))

    bin_dir = project_root / "bin"
    script_dir = Path(__file__).parent
    update_script = script_dir / "update_param.py"

    # --- Mode Dispatch ---
    
    if mode == "hdprm":
        # Run BHT_HDPRM, BH2_HDPRM, HTOF_HDPRM, BAC_HDPRM, KVC_HDPRM, T1_HDPRM, CVC_HDPRM, SAC3_HDPRM, SFV_HDPRM
        # Note: update_param.py supports these 9 detectors for 'hdprm'
        detectors = ["BHT", "BH2", "HTOF", "BAC", "KVC", "T1", "CVC", "SAC3", "SFV"]
        
        print(colored(f">>> Step 1: Running HDPRM analysis for {len(detectors)} detectors", "cyan"))
        for det in detectors:
            binary = bin_dir / f"{det}_HDPRM"
            if binary.exists():
                run_command(f"{binary} {input_root_file} {suffix}")
            else:
                print(colored(f"[Warning] Binary {binary.name} not found. Skipping.", "yellow"))
        
        print(colored(">>> Step 2: Updating Parameters (HDPRM)", "cyan"))
        run_command(f"python3 {update_script} {run_num} {suffix} hdprm")
        
    elif mode == "t0":
        # Run T0_Offset -> update_param.py t0
        executable = bin_dir / "T0_Offset"
        if not executable.exists():
            print(colored(f"[Error] T0_Offset not found. Please compile.", "red"))
            sys.exit(1)
            
        print(colored(">>> Step 1: Running T0_Offset", "cyan"))
        run_command(f"{executable} {input_root_file} {suffix}")
        
        print(colored(">>> Step 2: Updating Parameters (T0)", "cyan"))
        run_command(f"python3 {update_script} {run_num} {suffix} t0")
        
    elif mode == "hdphc":
        # Run BHT_PHC, BH2_PHC, HTOF_PHC, T1_PHC, CVC_PHC
        # Note: update_param.py supports these 5 detectors for 'hdphc'
        detectors = ["BHT", "BH2", "HTOF", "T1", "CVC"]
        
        print(colored(f">>> Step 1: Running PHC analysis for {len(detectors)} detectors", "cyan"))
        for det in detectors:
            binary = bin_dir / f"{det}_PHC"
            if binary.exists():
                run_command(f"{binary} {input_root_file} {suffix}")
            else:
                print(colored(f"[Warning] Binary {binary.name} not found. Skipping.", "yellow"))
        
        print(colored(">>> Step 2: Updating Parameters (HDPHC)", "cyan"))
        run_command(f"python3 {update_script} {run_num} {suffix} hdphc")

    print(colored(f"\n[DONE] Hodo Calibration Complete for mode: {mode}", "green", attrs=["bold"]))

if __name__ == "__main__":
    main()
