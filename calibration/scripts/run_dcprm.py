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
        description=(
            "Automate BC (DC) Calibration for T0, Drift (ST), Residual, and Resolution. "
            "reso: Pull Gauss width p -> Res *= p "
            "(iterate: rebuild Bc DST with #define FillPullExclusive 1 -> reso -> DST). "
            "Use --debug to evaluate pull without updating Res."
        )
    )
    parser.add_argument("run_num", type=int, help="Run Number")
    parser.add_argument(
        "mode",
        type=str,
        choices=["t0", "drift", "resi", "reso"],
        help="Tuning Mode (reso: Pull p -> DCGEO Res *= p; iterate with DST)",
    )
    
    # Mutually exclusive group for detector selection
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--bcout', action="store_true", help='Set mode to BcOut calibration (default if no flag is set)')
    group.add_argument('--bcin', action="store_true", help='Set mode to BcIn calibration')
    group.add_argument(
        '--all',
        action="store_true",
        help='reso only: run BLC1+BLC2 Pull, average p, write same Res scale to both',
    )
    parser.add_argument('--kaon', action="store_true", help='Use Kaon (K) suffix instead of Pion (Pi)')
    parser.add_argument('--guess', action="store_true",
                        help='Use full-wire projection as Erfc fit seed for TDC (legacy; t0 mode only)')
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Skip parameter update; run analysis and produce PDFs only",
    )
    
    args = parser.parse_args()

    run_num = args.run_num
    mode = args.mode

    if args.all and mode != "reso":
        print(colored("[Error] --all is only supported for mode=reso", "red"))
        sys.exit(1)
    
    # Determined detector based on flags
    # Default is BcOut
    detector = "BcOut"
    if args.bcin:
        detector = "BcIn"
        
    suffix = "K" if args.kaon else "Pi"

    def resolve_input_root(det_name):
        root_filename = f"run{run_num:05d}_{det_name}.root"
        path = config.DATA_DIR / root_filename
        if not path.exists() and det_name == "BcOut":
            fallback = config.DATA_DIR / f"run{run_num:05d}_DC.root"
            if fallback.exists():
                print(colored(f"[INFO] _BcOut.root not found, using _DC.root instead.", "yellow"))
                return fallback
        return path

    bin_dir = project_root / "bin"
    script_dir = Path(__file__).parent
    update_script = script_dir / "update_param.py"

    if args.all:
        # reso --all: both chambers
        input_roots = []
        for det_name in ("BcIn", "BcOut"):
            path = resolve_input_root(det_name)
            if not path.exists():
                print(colored(f"[Error] Symlink/File not found: {path}", "red"))
                print(f"Ensure create_runlist.py was run for {det_name}.")
                sys.exit(1)
            input_roots.append((det_name, path))
        print(colored("[INFO] Mode: reso --all (Suffix: {})".format(suffix), "green"))
        for det_name, path in input_roots:
            print(colored(f"[INFO] Input {det_name}: {path}", "green"))
        if args.debug:
            print(colored("[INFO] DEBUG mode: parameter update will be skipped", "yellow"))
    else:
        # 1. Locate the Input Root File (Symbolic link in DATA_DIR)
        input_root_file = resolve_input_root(detector)
             
        if not input_root_file.exists():
             print(colored(f"[Error] Symlink/File not found: {input_root_file}", "red"))
             print(f"Ensure create_runlist.py was run for {detector}.")
             sys.exit(1)
             
        print(colored(f"[INFO] Using Input File: {input_root_file}", "green"))
        print(colored(f"[INFO] Target Detector: {detector}", "green"))
        print(colored(f"[INFO] Mode: {mode} (Suffix: {suffix})", "green"))
        if args.debug:
            print(colored("[INFO] DEBUG mode: parameter update will be skipped", "yellow"))

    # --- Mode Dispatch ---
    
    if mode == "t0":
        # Run BLC_TDC -> udpate_param.py dctdc
        executable = bin_dir / "BLC_TDC"
        if not executable.exists():
            print(colored(f"[Error] BLC_TDC not found. Please compile.", "red"))
            sys.exit(1)
            
        print(colored(">>> Step 1: Running BLC_TDC", "cyan"))
        guess_flag = " --guess" if args.guess else ""
        run_command(f"{executable} {input_root_file} {suffix}{guess_flag}")
        
        if not args.debug:
            print(colored(">>> Step 2: Updating Parameters (DCTDC)", "cyan"))
            det_flag = "--bcin" if args.bcin else "--bcout"
            run_command(f"python3 {update_script} {run_num} {suffix} dctdc {det_flag}")
        else:
            print(colored(">>> [DEBUG] Skipping parameter update", "yellow"))
        
    elif mode == "drift":
        # Run BLC_DRFT -> (No update_param needed, it updates ROOT param file directly)
        executable = bin_dir / "BLC_DRFT"
        if not executable.exists():
            print(colored(f"[Error] BLC_DRFT not found. Please compile.", "red"))
            sys.exit(1)
            
        print(colored(">>> Step 1: Running BLC_DRFT", "cyan"))
        debug_flag = " --debug" if args.debug else ""
        run_command(f"{executable} {input_root_file} {suffix}{debug_flag}")
        
        if not args.debug:
            print(colored(">>> (Drift parameters are updated directly in param/DCDRFT/)", "yellow"))
        else:
            print(colored(">>> [DEBUG] Skipping parameter update", "yellow"))
        
    elif mode == "resi":
        # Run BLC_residual -> update_param.py residual
        executable = bin_dir / "BLC_residual"
        if not executable.exists():
            print(colored(f"[Error] BLC_residual not found. Please compile.", "red"))
            sys.exit(1)
            
        print(colored(">>> Step 1: Running BLC_residual", "cyan"))
        run_command(f"{executable} {input_root_file} {suffix}")
        
        if not args.debug:
            print(colored(">>> Step 2: Updating Parameters (Residual -> DCGEO Ofs)", "cyan"))
            # Using 'residual' type which has our special additive logic
            det_flag = "--bcin" if args.bcin else "--bcout"
            run_command(f"python3 {update_script} {run_num} {suffix} residual {det_flag}")
        else:
            print(colored(">>> [DEBUG] Skipping parameter update", "yellow"))

    elif mode == "reso":
        # BLC_pull: Pull Gauss width p; update DCGEO Res *= p (iterate with DST)
        executable = bin_dir / "BLC_pull"
        if not executable.exists():
            print(colored(f"[Error] BLC_pull not found. Please compile.", "red"))
            sys.exit(1)

        if args.all:
            print(colored(">>> Step 1: Running BLC_pull (BcIn + BcOut)", "cyan"))
            for det_name, path in input_roots:
                print(colored(f">>>   {det_name}", "cyan"))
                run_command(f"{executable} {path} {suffix}")
            if not args.debug:
                print(colored(
                    ">>> Step 2: Updating Parameters (Pull ALL -> DCGEO Res *= p)",
                    "cyan",
                ))
                run_command(f"python3 {update_script} {run_num} {suffix} resolution --all")
            else:
                print(colored(">>> [DEBUG] Skipping parameter update", "yellow"))
        else:
            print(colored(">>> Step 1: Running BLC_pull", "cyan"))
            run_command(f"{executable} {input_root_file} {suffix}")

            if not args.debug:
                print(colored(">>> Step 2: Updating Parameters (Pull -> DCGEO Res *= p)", "cyan"))
                det_flag = "--bcin" if args.bcin else "--bcout"
                run_command(f"python3 {update_script} {run_num} {suffix} resolution {det_flag}")
            else:
                print(colored(">>> [DEBUG] Skipping parameter update", "yellow"))

    print(colored(f"\n[DONE] Tuning Complete for mode: {mode}", "green", attrs=["bold"]))

if __name__ == "__main__":
    main()
