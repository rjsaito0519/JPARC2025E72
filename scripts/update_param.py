#!/usr/bin/env python3

# argparse.ArgumentParserクラスをインスタンス化して、説明等を引数として渡す
# ---------------------------------------------------------------------------
import argparse
parser = argparse.ArgumentParser(
    prog="update_param",
    usage="python3 update_param.py <run_num> <suffix> <param_type> <--ftof>", # プログラムの利用方法
    description="", # ヘルプの前に表示
    epilog="end", # ヘルプの後に表示
    add_help=True, # -h/–-helpオプションの追加
)
parser.add_argument("run_num", type=int, help="Input run number")
parser.add_argument("suffix", type=str, help="Input suffix")
parser.add_argument("param_type", type=str, help="Input parameter type")
parser.add_argument('--ftof', action="store_true", help='Set check to True')
args = parser.parse_args()
# ---------------------------------------------------------------------------

import sys
if args.suffix not in ["K", "Pi"]:
    print("suffix should be K or Pi")
    sys.exit()

if args.param_type not in ["hdprm", "t0", "hdphc", "dctdc", "residual"]:
    print("param_type should be hdprm or t0 or hdphc or dctdc or residual")
    sys.exit()

import os
import shutil
import numpy as np
from termcolor import colored

def report_status(do_succeeded, counter_name):
    if do_succeeded:
        print(colored(f"[OK]   {counter_name} successfully updated", "green"))
    else:
        print(colored(f"[FAIL] {counter_name} update failed", "red"))

import update_hdprm
import update_phc
import phc_conf
import hdprm_conf
import update_dctdc
import update_residual

# prepare param file
# ---------------------------------------------------------------------------
import conf

SUB_DIR = "e72"

hdprm_dir = f"{conf.param_dir}/HDPRM"
hdprm_target_file = f"{hdprm_dir}/{SUB_DIR}/HodoParam_run{args.run_num:0=5}_{args.suffix}"
if not os.path.isfile(hdprm_target_file):
    shutil.copy(f"{hdprm_dir}/HodoParam_0", hdprm_target_file)

hdphc_dir = f"{conf.param_dir}/HDPHC"
hdphc_target_file = f"{hdphc_dir}/{SUB_DIR}/HodoPHCParam_run{args.run_num:0=5}_{args.suffix}"
if not os.path.isfile(hdphc_target_file):
    shutil.copy(f"{hdphc_dir}/HodoPHCParam_0", hdphc_target_file)

dctdc_dir = f"{conf.param_dir}/DCTDC"
dctdc_target_file = f"{dctdc_dir}/{SUB_DIR}/DCTdcParam_run{args.run_num:0=5}_{args.suffix}"
if not os.path.isfile(dctdc_target_file):
    shutil.copy(f"{dctdc_dir}/DCTdcParam_0", dctdc_target_file)

dcgeo_dir = f"{conf.param_dir}/DCGEO"
dcgeo_target_file = f"{dcgeo_dir}/{SUB_DIR}/DCGeomParam_run{args.run_num:0=5}_{args.suffix}"
if not os.path.isfile(dcgeo_target_file):
    shutil.copy(f"{dcgeo_dir}/DCGeomParam_e73_2024_0", dcgeo_target_file)

# update param file
# ---------------------------------------------------------------------------
if args.param_type == "hdprm":
    detectors = ["BHT", "T0"]
    for det in detectors:
        root_file = os.path.join(
            conf.data_dir,
            f"root/run{args.run_num:0=5}_{det}_HDPRM_{args.suffix}.root"
        )
        if not os.path.exists(root_file):
            print(colored(f"[SKIP] {det}: file not found", "yellow"))
            continue

        limits = [-np.inf, np.inf]
        if f"{args.run_num:0=5}_{args.suffix}_{det.lower()}" in hdprm_conf.limits_dict.keys():
            limits = hdprm_conf.limits_dict[f"{args.run_num:0=5}_{args.suffix}_{det.lower()}"]
        data = update_hdprm.make_dictdata(root_file, limits)
        do_succeeded = update_hdprm.update_file(hdprm_target_file, data)
        report_status(do_succeeded, det)

elif args.param_type == "t0":
    # -- T0 -----
    data = update_hdprm.make_dictdata(os.path.join(conf.data_dir, f"root/run{args.run_num:0=5}_T0_Offset_{args.suffix}.root"), [-np.inf, np.inf], is_t0_offset = True)
    do_succeeded = update_hdprm.update_file(hdprm_target_file, data)
    report_status(do_succeeded, "T0")

elif args.param_type == "hdphc":
    detectors = ["BHT", "T0"]
    for det in detectors:
        root_file = os.path.join(
            conf.data_dir,
            f"root/run{args.run_num:0=5}_{det}_PHC_{args.suffix}.root"
        )
        if not os.path.exists(root_file):
            print(colored(f"[SKIP] {det}: file not found", "yellow"))
            continue

        limits = [-np.inf, np.inf]
        if f"{args.run_num:0=5}_{args.suffix}_{det.lower()}" in phc_conf.limits_dict.keys():
            limits = phc_conf.limits_dict[f"{args.run_num:0=5}_{args.suffix}_{det.lower()}"]
        data = update_phc.make_dictdata(root_file, limits)
        do_succeeded = update_phc.update_file(hdphc_target_file, data)
        report_status(do_succeeded, det)

elif args.param_type == "dctdc":
    # # -- BLC1 -----
    # data = update_dctdc.make_dictdata(os.path.join(conf.data_dir, f"root/run{args.run_num:0=5}_BLC1_TDC_{args.suffix}.root"))
    # do_succeeded = update_dctdc.update_file(dctdc_target_file, data)
    # report_status(do_succeeded, "BLC1")
    
    # -- BLC2 -----
    data = update_dctdc.make_dictdata(os.path.join(conf.data_dir, f"root/run{args.run_num:0=5}_BLC2_TDC_{args.suffix}.root"))
    do_succeeded = update_dctdc.update_file(dctdc_target_file, data)
    report_status(do_succeeded, "BLC2")

elif args.param_type == "residual":
    # # -- BLC1 -----
    # data = update_dctdc.make_dictdata(os.path.join(conf.data_dir, f"root/run{args.run_num:0=5}_BLC1_TDC_{args.suffix}.root"))
    # do_succeeded = update_dctdc.update_file(dctdc_target_file, data)
    # report_status(do_succeeded, "BLC1")
    
    # -- BLC2 -----
    data = update_residual.make_dictdata(os.path.join(conf.data_dir, f"root/run{args.run_num:0=5}_BLC2_residual_{args.suffix}.root"))
    do_succeeded = update_residual.update_file(dcgeo_target_file, data)
    report_status(do_succeeded, "BLC2")

