#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path
from termcolor import colored

# Add shared library path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config

def run_command(cmd):
    print(colored(f"Executing: {cmd}", "cyan"))
    try:
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end="")
        process.wait()
        return process.returncode == 0
    except Exception as e:
        print(colored(f"Error: {e}", "red"))
        return False

def main():
    parser = argparse.ArgumentParser(
        prog="run_hodo",
        description="Automate Hodoscope calibration (HDPRM, T0, HDPHC) for T110 Branch.",
        add_help=True,
    )
    parser.add_argument("mode", choices=["hdprm", "t0", "hdphc"], help="Calibration mode")
    parser.add_argument("input", type=str, help="Input ROOT file or Run number")
    parser.add_argument("suffix", choices=["K", "Pi"], help="Particle suffix (K or Pi)")
    parser.add_argument("--kaon", action="store_true", help="Shortcut for 'K' suffix (Particle name determined by this)")
    
    args = parser.parse_args()
    suffix = "K" if args.kaon else args.suffix
    
    # 1. Determine run number and ROOT file
    input_path = Path(args.input)
    run_num = None
    input_root_file = None
    
    if input_path.exists() and input_path.is_file():
        input_root_file = input_path.resolve()
        import re
        match = re.search(r"run(\d{5})", str(input_root_file))
        if match: run_num = int(match.group(1))
    elif args.input.isdigit():
        run_num = int(args.input)
        # Search priority: Decode -> Scratch -> Output
        search_paths = [
            config.DECODE_DIR / f"run{run_num:05d}" / f"run{run_num:05d}_Hodo.root",
            config.SCRATCH_DIR / f"run{run_num:05d}" / f"run{run_num:05d}_Hodo.root",
            config.OUTPUT_DIR / "root" / f"run{run_num:05d}" / f"run{run_num:05d}_Hodo.root"
        ]
        for p in search_paths:
            if p.exists():
                input_root_file = p
                break
    
    if not input_root_file or not input_root_file.exists():
        print(colored(f"[Error] Could not find input ROOT file for Run {run_num}", "red"))
        sys.exit(1)

    print(colored(f"\n>>> Starting Hodo Calibration ({args.mode.upper()}) for Run {run_num} ({suffix}) <<<", "green", attrs=["bold"]))
    print(f"Input: {input_root_file}\n")

    bin_dir = project_root / "bin"
    update_script = project_root / "calibration" / "scripts" / "update_param.py"
    mode = args.mode

    # --- Mode Dispatch ---
    if mode == "hdprm":
        # T110 standard detectors for HDPRM
        detectors = ["BHT", "T0", "BH2", "BAC", "KVC", "T1", "CVC", "SAC3", "SFV"]
        
        print(colored(f">>> Step 1: Running HDPRM analysis", "cyan"))
        for det in detectors:
            binary = bin_dir / f"{det}_HDPRM"
            if binary.exists():
                run_command(f"{binary} {input_root_file} {suffix}")
            else:
                pass # Skip silently for robustness
        
        print(colored(">>> Step 2: Updating Parameters (HDPRM)", "cyan"))
        run_command(f"python3 {update_script} {run_num} {suffix} hdprm")
        
    elif mode == "t0":
        executable = bin_dir / "T0_Offset"
        if not executable.exists():
            print(colored(f"[Error] T0_Offset binary not found.", "red"))
            sys.exit(1)
            
        print(colored(">>> Step 1: Running T0_Offset", "cyan"))
        run_command(f"{executable} {input_root_file} {suffix}")
        
        print(colored(">>> Step 2: Updating Parameters (T0)", "cyan"))
        run_command(f"python3 {update_script} {run_num} {suffix} t0")
        
    elif mode == "hdphc":
        # T110 standard detectors for HDPHC
        detectors = ["BHT", "T0", "BH2", "T1", "CVC"]
        
        print(colored(f">>> Step 1: Running PHC analysis", "cyan"))
        for det in detectors:
            binary = bin_dir / f"{det}_PHC"
            if binary.exists():
                run_command(f"{binary} {input_root_file} {suffix}")
        
        print(colored(">>> Step 2: Updating Parameters (HDPHC)", "cyan"))
        run_command(f"python3 {update_script} {run_num} {suffix} hdphc")

    print(colored(f"\n[DONE] Hodo {mode.upper()} calibration finished for Run {run_num}.", "green", attrs=["bold"]))

if __name__ == "__main__":
    main()
