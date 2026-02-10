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
        description="Automate BC (DC) Calibration for T0, Drift (ST), and Residual."
    )
    parser.add_argument("run_num", type=int, help="Run Number")
    parser.add_argument("mode", type=str, choices=["t0", "drift", "resi"], help="Tuning Mode")
    
    # Mutually exclusive group for detector selection
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--bcout', action="store_true", help='Set mode to BcOut calibration (default if no flag is set)')
    group.add_argument('--bcin', action="store_true", help='Set mode to BcIn calibration')
    parser.add_argument('--kaon', action="store_true", help='Use Kaon (K) suffix instead of Pion (Pi)')
    
    args = parser.parse_args()

    run_num = args.run_num
    mode = args.mode
    
    # Determined detector based on flags
    # Default is BcOut
    detector = "BcOut"
    if args.bcin:
        detector = "BcIn"
        
    suffix = "K" if args.kaon else "Pi"

    # 1. Locate the Input Root File (Symbolic link in DATA_DIR)
    # create_runlist.py creates runXXXXX_{detector}.root
    
    root_filename = f"run{run_num:05d}_{detector}.root"
    input_root_file = config.DATA_DIR / root_filename
    
    # Fallback to _DC if _BcOut not found (compatibility)
    if not input_root_file.exists():
         # Only fallback to DC for BcOut as that was the old naming
         if detector == "BcOut":
             fallback = config.DATA_DIR / f"run{run_num:05d}_DC.root"
             if fallback.exists():
                 print(colored(f"[INFO] _BcOut.root not found, using _DC.root instead.", "yellow"))
                 input_root_file = fallback
         
    if not input_root_file.exists():
         print(colored(f"[Error] Symlink/File not found: {input_root_file}", "red"))
         print(f"Ensure create_runlist.py was run for {detector}.")
         sys.exit(1)
         
    print(colored(f"[INFO] Using Input File: {input_root_file}", "green"))
    print(colored(f"[INFO] Target Detector: {detector}", "green"))
    print(colored(f"[INFO] Mode: {mode} (Suffix: {suffix})", "green"))

    bin_dir = project_root / "bin"
    script_dir = Path(__file__).parent
    update_script = script_dir / "update_param.py"

    # --- Mode Dispatch ---
    
    if mode == "t0":
        # Run BLC_TDC -> udpate_param.py dctdc
        executable = bin_dir / "BLC_TDC"
        if not executable.exists():
            print(colored(f"[Error] BLC_TDC not found. Please compile.", "red"))
            sys.exit(1)
            
        print(colored(">>> Step 1: Running BLC_TDC", "cyan"))
        run_command(f"{executable} {input_root_file} {suffix}")
        
        print(colored(">>> Step 2: Updating Parameters (DCTDC)", "cyan"))
        det_flag = "--bcin" if args.bcin else "--bcout"
        run_command(f"python3 {update_script} {run_num} {suffix} dctdc {det_flag}")
        
    elif mode == "drift":
        # Run BLC_DRFT -> (No update_param needed, it updates ROOT param file directly)
        executable = bin_dir / "BLC_DRFT"
        if not executable.exists():
            print(colored(f"[Error] BLC_DRFT not found. Please compile.", "red"))
            sys.exit(1)
            
        print(colored(">>> Step 1: Running BLC_DRFT", "cyan"))
        run_command(f"{executable} {input_root_file} {suffix}")
        
        print(colored(">>> (Drift parameters are updated directly in param/DCDRFT/)", "yellow"))
        
    elif mode == "resi":
        # Run BLC_residual -> update_param.py residual
        executable = bin_dir / "BLC_residual"
        if not executable.exists():
            print(colored(f"[Error] BLC_residual not found. Please compile.", "red"))
            sys.exit(1)
            
        print(colored(">>> Step 1: Running BLC_residual", "cyan"))
        run_command(f"{executable} {input_root_file} {suffix}")
        
        print(colored(">>> Step 2: Updating Parameters (Residual -> DCGEO)", "cyan"))
        # Using 'residual' type which has our special additive logic
        det_flag = "--bcin" if args.bcin else "--bcout"
        run_command(f"python3 {update_script} {run_num} {suffix} residual {det_flag}")

    print(colored(f"\n[DONE] Tuning Complete for mode: {mode}", "green", attrs=["bold"]))

if __name__ == "__main__":
    main()
