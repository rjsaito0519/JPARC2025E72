# argparse.ArgumentParserクラスをインスタンス化して、説明等を引数として渡す
# ---------------------------------------------------------------------------
import argparse
parser = argparse.ArgumentParser(
    prog="update_param",
    usage="python3 update_param.py <run_num> <suffix> <param_type>", # プログラムの利用方法
    description="", # ヘルプの前に表示
    epilog="end", # ヘルプの後に表示
    add_help=True, # -h/–-helpオプションの追加
)
parser.add_argument("run_num", type=int, help="Input run number")
parser.add_argument("suffix", type=str, help="Input suffix")
parser.add_argument("param_type", type=str, help="Input parameter type")
args = parser.parse_args()
# ---------------------------------------------------------------------------

import sys
if args.suffix not in ["K", "Pi"]:
    print("suffix should be K or Pi")
    sys.exit()

if args.param_type not in ["hdprm", "t0", "hdphc", "dctdc" "residual"]:
    print("param_type should be hdprm or t0 or hdphc or dctdc or residual")
    sys.exit()

import os
import shutil
import numpy as np

def report_status(do_succeeded, counter_name):
    if do_succeeded:
        print(f"\033[92m✅ {counter_name} {args.param_type} successfully updated\033[0m")
    else:
        print(f"\033[91m❌ Failed to update {counter_name} {args.param_type}\033[0m")

import update_hdprm
import update_phc
import phc_conf
import update_dctdc
import update_residual
script_dir = os.path.dirname(os.path.abspath(__file__))

# prepare param file
# ---------------------------------------------------------------------------
import conf

SUB_DIR = "e72"

hdprm_dir = f"{conf.param_dir}/HDPRM"
hdprm_target_file = f"{hdprm_dir}/{SUB_DIR}/HodoParam_run{args.run_num:0=5}_{args.suffix}"
if not os.path.isfile(hdprm_target_file):
    shutil.copy(f"{hdprm_dir}/HodoParam_1", hdprm_target_file)

hdphc_dir = f"{conf.param_dir}/HDPHC"
hdphc_target_file = f"{hdphc_dir}/{SUB_DIR}/HodoPHCParam_run{args.run_num:0=5}_{args.suffix}"
if not os.path.isfile(hdphc_target_file):
    shutil.copy(f"{hdphc_dir}/HodoPHCParam_1", hdphc_target_file)

dctdc_dir = f"{conf.param_dir}/DCTDC"
dctdc_target_file = f"{dctdc_dir}/{SUB_DIR}/DCTdcParam_run{args.run_num:0=5}_{args.suffix}"
if not os.path.isfile(dctdc_target_file):
    shutil.copy(f"{dctdc_dir}/DCTdcParam_1", dctdc_target_file)

dcgeo_dir = f"{conf.param_dir}/DCGEO"
dcgeo_target_file = f"{dcgeo_dir}/{SUB_DIR}/DCGeomParam_run{args.run_num:0=5}_{args.suffix}"
if not os.path.isfile(dcgeo_target_file):
    shutil.copy(f"{dcgeo_dir}/DCGeomParam_e72_20251120", dcgeo_target_file)

# update param file
# ---------------------------------------------------------------------------
if args.param_type == "hdprm":
    # -- BHT -----
    data = update_hdprm.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BHT_HDPRM_{args.suffix}.root"))
    do_succeeded = update_hdprm.update_file(hdprm_target_file, data)
    report_status(do_succeeded, "BHT")

    # -- BH2 -----
    data = update_hdprm.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BH2_HDPRM_{args.suffix}.root"))
    do_succeeded = update_hdprm.update_file(hdprm_target_file, data)
    report_status(do_succeeded, "BH2")

    # -- HTOF -----
    data = update_hdprm.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_HTOF_HDPRM_{args.suffix}.root"))
    do_succeeded = update_hdprm.update_file(hdprm_target_file, data)
    report_status(do_succeeded, "HTOF")

elif args.param_type == "t0":
    # -- T0 -----
    data = update_hdprm.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_T0_Offset_{args.suffix}.root"), is_t0_offset = True)
    do_succeeded = update_hdprm.update_file(hdprm_target_file, data)
    report_status(do_succeeded, "BH2")

elif args.param_type == "hdphc":
    # -- BHT -----
    limits = [-np.inf, np.inf]
    if f"{args.run_num:0=5}_{args.suffix}_bht" in phc_conf.limits_dict.keys():
        limits = phc_conf.limits_dict[f"{args.run_num:0=5}_{args.suffix}_bht"]
    data = update_phc.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BHT_PHC_{args.suffix}.root"), limits)
    do_succeeded = update_phc.update_file(hdphc_target_file, data)
    report_status(do_succeeded, "BHT")

    # -- BH2 -----
    limits = [-np.inf, np.inf]
    if f"{args.run_num:0=5}_{args.suffix}_bh2" in phc_conf.limits_dict.keys():
        limits = phc_conf.limits_dict[f"{args.run_num:0=5}_{args.suffix}_bh2"]
    data = update_phc.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BH2_PHC_{args.suffix}.root"), limits)
    do_succeeded = update_phc.update_file(hdphc_target_file, data)
    report_status(do_succeeded, "BH2")

    # # -- HTOF -----
    # limits = [-np.inf, np.inf]
    # if f"{args.run_num:0=5}_{args.suffix}_htof" in phc_conf.limits_dict.keys():
    #     limits = phc_conf.limits_dict[f"{args.run_num:0=5}_{args.suffix}_htof"]
    # data = update_phc.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_HTOF_PHC_{args.suffix}.root"), limits)
    # do_succeeded = update_phc.update_file(hdphc_target_file, data)
    # report_status(do_succeeded, "HTOF")


elif args.param_type == "dctdc":
    # -- BLC1 -----
    data = update_dctdc.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BLC1_TDC_{args.suffix}.root"))
    do_succeeded = update_dctdc.update_file(dctdc_target_file, data)
    report_status(do_succeeded, "BLC1")
    
    # -- BLC2 -----
    data = update_dctdc.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BLC2_TDC_{args.suffix}.root"))
    do_succeeded = update_dctdc.update_file(dctdc_target_file, data)
    report_status(do_succeeded, "BLC2")

elif args.param_type == "residual":
    # # -- BLC1 -----
    # data = update_dctdc.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BLC1_TDC_{args.suffix}.root"))
    # do_succeeded = update_dctdc.update_file(dctdc_target_file, data)
    # report_status(do_succeeded, "BLC1")
    
    # -- BLC2 -----
    data = update_residual.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BLC2_residual_{args.suffix}.root"))
    do_succeeded = update_residual.update_file(dcgeo_target_file, data)
    report_status(do_succeeded, "BLC2")

