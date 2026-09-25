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
        description="Automate Hodoscope Calibration for HDPRM, T0 Offset, HDPHC, and HTOF (TDC/ADC/PHC)."
    )
    parser.add_argument("run_nums", type=int, nargs="+",
                        help="Run Number(s). Several runs are merged (htof only); "
                             "results/params are labeled with the first run")
    parser.add_argument("mode", type=str, choices=["hdprm", "htof", "t0", "hdphc"], help="Calibration Mode")
    parser.add_argument('--kaon', action="store_true", help='Use Kaon (K) suffix instead of Pion (Pi)')
    parser.add_argument('--ftof', action="store_true", help='Include FTOF-related detectors (CVC, SFV, SAC3 for hdprm; CVC for hdphc)')
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Skip parameter update; run analysis and produce PDFs only",
    )

    args = parser.parse_args()

    run_nums = args.run_nums
    run_num = run_nums[0]  # representative run (output/param labeling)
    mode = args.mode
    suffix = "K" if args.kaon else "Pi"

    if len(run_nums) > 1 and mode != "htof":
        print(colored(f"[Error] Multiple runs are supported only for htof (got mode={mode}).", "red"))
        sys.exit(1)

    if mode == "htof" and args.kaon:
        print(colored("[Error] htof only supports Pi (TPC dE/dx pion bit set, proton bit excluded); "
                       "--kaon is not implemented.", "red"))
        sys.exit(1)

    # 1. Locate the Input Root File(s) (Symbolic link in DATA_DIR)
    # create_runlist.py creates runXXXXX_Hodo.root
    hodo_root_files = []
    for r in run_nums:
        hodo_root_file = config.DATA_DIR / f"run{r:05d}_Hodo.root"
        if not hodo_root_file.exists():
            print(colored(f"[Error] Symlink/File not found: {hodo_root_file}", "red"))
            print(f"Ensure create_runlist.py was run for Hodo.")
            sys.exit(1)
        hodo_root_files.append(hodo_root_file)
    input_root_file = hodo_root_files[0]

    htof_dst_files = []
    if mode == "htof":
        # HTOF_Calib also needs the DstTPCHelixHTOF output (dE/dx-tagged TPC helix x
        # HTOF match) for the pion-selected MIP histogram; Hodo.root alone supplies
        # the all-event pedestal histograms. Not yet wired into
        # create_runlist.py/dst_create_runlist.py; produce it manually for now
        # (bin/DstTPCHelixHTOF <conf> <TPCHelix.root> <Hodo.root> run{run}_HTOF.root).
        for r in run_nums:
            htof_dst_file = config.DATA_DIR / f"run{r:05d}_HTOF.root"
            if not htof_dst_file.exists():
                print(colored(f"[Error] Symlink/File not found: {htof_dst_file}", "red"))
                print("Run bin/DstTPCHelixHTOF for this run first (dst_create_runlist.py "
                      "integration is not yet available).")
                sys.exit(1)
            htof_dst_files.append(htof_dst_file)
            print(colored(f"[INFO] Using DstTPCHelixHTOF File: {htof_dst_file}", "green"))

    for f in hodo_root_files:
        print(colored(f"[INFO] Using Input File: {f}", "green"))
    if len(run_nums) > 1:
        print(colored(f"[INFO] Merging {len(run_nums)} runs; results/params labeled with run {run_num:05d}", "green"))
    print(colored(f"[INFO] Mode: {mode} (Suffix: {suffix})", "green"))
    if args.ftof:
        print(colored("[INFO] FTOF detectors included", "green"))
    if args.debug:
        print(colored("[INFO] DEBUG mode: parameter update will be skipped", "yellow"))

    bin_dir = project_root / "bin"
    script_dir = Path(__file__).parent
    update_script = script_dir / "update_param.py"

    # --- Mode Dispatch ---
    
    if mode == "hdprm":
        # Run BHT_HDPRM, BH2_HDPRM, BAC_HDPRM, KVC_HDPRM, T1_HDPRM, ...
        # HTOF is handled by mode "htof" (HTOF_Calib).
        detectors = ["BHT", "BH2", "BAC", "KVC", "T1"] + (["CVC", "SFV", "SAC3"] if args.ftof else [])
        
        print(colored(f">>> Step 1: Running HDPRM analysis for {len(detectors)} detectors", "cyan"))
        for det in detectors:
            binary = bin_dir / f"{det}_HDPRM"
            if binary.exists():
                run_command(f"{binary} {input_root_file} {suffix}")
            else:
                print(colored(f"[Warning] Binary {binary.name} not found. Skipping.", "yellow"))
        
        if not args.debug:
            print(colored(">>> Step 2: Updating Parameters (HDPRM)", "cyan"))
            run_command(f"python3 {update_script} {run_num} {suffix} hdprm")
        else:
            print(colored(">>> [DEBUG] Skipping parameter update", "yellow"))
        
    elif mode == "htof":
        # Run HTOF_Calib (TDC + ADC pedestal/MIP + PHC in one pass; PHC uses DeltaE
        # recomputed from raw ADC with the new ADC parameters, so no re-decoding is
        # needed in between) -> update_param.py htofprm (HDPRM: TDC+ADC) and
        # htofphc (HDPHC). One PDF per segment: run{run}_HTOF_Calib_{suffix}.pdf.
        executable = bin_dir / "HTOF_Calib"
        if not executable.exists():
            print(colored(f"[Error] HTOF_Calib not found. Please compile.", "red"))
            sys.exit(1)

        print(colored(">>> Step 1: Running HTOF_Calib", "cyan"))
        pairs = " ".join(f"{h} {d}" for h, d in zip(hodo_root_files, htof_dst_files))
        run_command(f"{executable} {pairs} {suffix}")

        if not args.debug:
            print(colored(">>> Step 2: Updating Parameters (HTOF HDPRM: TDC + ADC)", "cyan"))
            run_command(f"python3 {update_script} {run_num} {suffix} htofprm")
            print(colored(">>> Step 3: Updating Parameters (HTOF HDPHC)", "cyan"))
            run_command(f"python3 {update_script} {run_num} {suffix} htofphc")
        else:
            print(colored(">>> [DEBUG] Skipping parameter update", "yellow"))

    elif mode == "t0":
        # Run T0_Offset -> update_param.py t0
        executable = bin_dir / "T0_Offset"
        if not executable.exists():
            print(colored(f"[Error] T0_Offset not found. Please compile.", "red"))
            sys.exit(1)
            
        print(colored(">>> Step 1: Running T0_Offset", "cyan"))
        run_command(f"{executable} {input_root_file} {suffix}")
        
        if not args.debug:
            print(colored(">>> Step 2: Updating Parameters (T0)", "cyan"))
            run_command(f"python3 {update_script} {run_num} {suffix} t0")
        else:
            print(colored(">>> [DEBUG] Skipping parameter update", "yellow"))
        
    elif mode == "hdphc":
        # Run BHT_PHC, BH2_PHC, T1_PHC, CVC_PHC (HTOF PHC is handled by mode "htof")
        # Note: update_param.py supports BHT, BH2, T1, CVC for 'hdphc'
        detectors = ["BHT", "BH2", "T1"] + (["CVC"] if args.ftof else [])
        
        print(colored(f">>> Step 1: Running PHC analysis for {len(detectors)} detectors", "cyan"))
        for det in detectors:
            binary = bin_dir / f"{det}_PHC"
            if binary.exists():
                run_command(f"{binary} {input_root_file} {suffix}")
            else:
                print(colored(f"[Warning] Binary {binary.name} not found. Skipping.", "yellow"))
        
        if not args.debug:
            print(colored(">>> Step 2: Updating Parameters (HDPHC)", "cyan"))
            run_command(f"python3 {update_script} {run_num} {suffix} hdphc")
        else:
            print(colored(">>> [DEBUG] Skipping parameter update", "yellow"))

    print(colored(f"\n[DONE] Hodo Calibration Complete for mode: {mode}", "green", attrs=["bold"]))

if __name__ == "__main__":
    main()
