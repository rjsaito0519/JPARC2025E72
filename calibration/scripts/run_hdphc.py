#!/usr/bin/env python3

import argparse
import subprocess
import sys
import re
from pathlib import Path
from termcolor import colored

# Add shared library path
# Current file: <RepoRoot>/calibration/scripts/run_hdphc.py
# parent: scripts, parentx2: calibration, parentx3: RepoRoot
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config

def run_command(cmd, shell=False):
    print(colored(f"Executing: {' '.join(cmd) if isinstance(cmd, list) else cmd}", "cyan"))
    try:
        process = subprocess.Popen(cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end="")
        process.wait()
        return process.returncode == 0
    except Exception as e:
        print(colored(f"Error: {e}", "red"))
        return False

def extract_run_num(path_str):
    match = re.search(r"run(\d{5})", path_str)
    if match:
        return int(match.group(1))
    return None

def main():
    parser = argparse.ArgumentParser(
        prog="run_hdphc",
        description="Run multiple HDPHC analysis binaries for a given ROOT file or run number.",
        add_help=True,
    )
    parser.add_argument("input", type=str, help="Input ROOT file path or Run number")
    parser.add_argument("suffix", type=str, choices=["K", "Pi"], help="Particle suffix (K or Pi)")
    parser.add_argument("--ftof", action="store_true", help="Include Forward TOF (HTOF, CVC) analysis")
    parser.add_argument("--update", action="store_true", help="Automatically call update_param.py after analysis")
    
    args = parser.parse_args()
    
    # 1. Determine input file and run number
    input_path = Path(args.input)
    if input_path.exists() and input_path.is_file():
        hodo_decoded_file = input_path.resolve()
        run_num = extract_run_num(str(hodo_decoded_file))
    elif args.input.isdigit():
        run_num = int(args.input)
        # 1. Try combined file in DECODE_DIR
        hodo_decoded_file = config.DECODE_DIR / f"run{run_num:05d}" / f"run{run_num:05d}_Hodo.root"
        if not hodo_decoded_file.exists():
            # 2. Try split file in DECODE_DIR
            hodo_decoded_file = config.DECODE_DIR / f"run{run_num:05d}" / f"run{run_num:05d}_Hodo_0.root"
        
        # 3. Try OUTPUT_DIR if not found
        if not hodo_decoded_file.exists():
            hodo_decoded_file = config.OUTPUT_DIR / "root" / f"run{run_num:05d}" / f"run{run_num:05d}_Hodo.root"
            if not hodo_decoded_file.exists():
                hodo_decoded_file = config.OUTPUT_DIR / "root" / f"run{run_num:05d}" / f"run{run_num:05d}_Hodo_0.root"

    else:
        print(colored(f"Error: Input '{args.input}' is neither an existing file nor a run number.", "red"))
        sys.exit(1)
    
    if not hodo_decoded_file.exists():
        print(colored(f"Error: Could not find decoded Hodo file: {hodo_decoded_file}", "red"))
        sys.exit(1)

    if run_num is None:
        print(colored("Error: Run number unknown (path doesn't match 'runXXXXX').", "red"))
        sys.exit(1)

    # 2. Define standard HDPHC binaries
    # Directory structure: <RepoRoot>/bin/
    bin_dir = project_root / "bin"
    
    standard_detectors = ["BHT", "T0"]
    ftof_detectors = ["CVC"]
    
    detectors_to_run = standard_detectors
    if args.ftof:
        detectors_to_run += ftof_detectors
        
    print(colored(f"\n>>> Starting HDPHC Analysis for Run {run_num} ({args.suffix}) <<<", "green", attrs=["bold"]))
    print(colored(f"Input file: {hodo_decoded_file}\n", "white"))

    failed = []
    for det in detectors_to_run:
        binary_name = f"{det}_PHC"
        binary_path = bin_dir / binary_name
        
        if not binary_path.exists():
            print(colored(f"Skip: {binary_name} not found in {bin_dir}", "yellow"))
            continue

        print(colored(f"--- Processing {det} ---", "magenta"))
        success = run_command([str(binary_path), str(hodo_decoded_file), args.suffix])
        
        if not success:
            failed.append(det)

    # Summary
    print("\n" + "="*50)
    if not failed:
        print(colored("All requested HDPHC analyses completed successfully!", "green", attrs=["bold"]))
    else:
        print(colored(f"Analyses failed for: {', '.join(failed)}", "red", attrs=["bold"]))
    print("="*50 + "\n")

    # 3. Optional Update
    if args.update and not failed:
        # Directory structure: <RepoRoot>/calibration/scripts/update_param.py
        update_script = project_root / "calibration" / "scripts" / "update_param.py"
        run_command([sys.executable, str(update_script), str(run_num), args.suffix, "hdphc"])

if __name__ == "__main__":
    main()
