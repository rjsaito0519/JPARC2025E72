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

if args.param_type not in ["hdprm", "t0", "hdphc"]:
    print("param_type should be hdprm or t0 or hdphc")
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
script_dir = os.path.dirname(os.path.abspath(__file__))

# prepare param file
# ---------------------------------------------------------------------------
import conf

hdprm_dir = f"{conf.param_dir}/HDPRM/hodo"
hdprm_target_file = f"{hdprm_dir}/HodoParam_run{args.run_num:0=5}_{args.suffix}"
if not os.path.isfile(hdprm_target_file):
    shutil.copy(f"{hdprm_dir}/HodoParam_0", hdprm_target_file)

hdphc_dir = f"{conf.param_dir}/HDPHC/hodo"
hdphc_target_file = f"{hdphc_dir}/HodoPHCParam_run{args.run_num:0=5}_{args.suffix}"
if not os.path.isfile(hdphc_target_file):
    shutil.copy(f"{hdphc_dir}/HodoPHCParam_0", hdphc_target_file)

# update param file
# ---------------------------------------------------------------------------
if args.param_type == "hdprm":
    # -- BHT -----
    data = update_hdprm.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BHT_HDPRM_{args.suffix}.root"))
    do_succeeded = update_hdprm.update_file(hdprm_target_file, data)
    report_status(do_succeeded, "BHT")

    # # -- T0 -----
    # data = update_hdprm.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_T0_HDPRM_{args.suffix}.root"))
    # do_succeeded = update_hdprm.update_file(hdprm_target_file, data)
    # report_status(do_succeeded, "T0")

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
    report_status(do_succeeded, "T0")

elif args.param_type == "hdphc":
    # -- BHT -----
    limits = [-np.inf, np.inf]
    if f"{args.run_num:0=5}_{args.suffix}_bht" in phc_conf.limits_dict.keys():
        limits = phc_conf.limits_dict[f"{args.run_num:0=5}_{args.suffix}_bht"]
    data = update_phc.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BHT_PHC_{args.suffix}.root"), limits)
    do_succeeded = update_phc.update_file(hdphc_target_file, data)
    report_status(do_succeeded, "BHT")

    # # -- T0 -----
    # limits = [-np.inf, np.inf]
    # if f"{args.run_num:0=5}_{args.suffix}_t0" in phc_conf.limits_dict.keys():
    #     limits = phc_conf.limits_dict[f"{args.run_num:0=5}_{args.suffix}_t0"]
    # data = update_phc.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_T0_PHC_{args.suffix}.root"), limits)
    # do_succeeded = update_phc.update_file(hdphc_target_file, data)
    # report_status(do_succeeded, "T0")

    # -- BH2 -----
    limits = [-np.inf, np.inf]
    if f"{args.run_num:0=5}_{args.suffix}_bh2" in phc_conf.limits_dict.keys():
        limits = phc_conf.limits_dict[f"{args.run_num:0=5}_{args.suffix}_bh2"]
    data = update_phc.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_BH2_PHC_{args.suffix}.root"), limits)
    do_succeeded = update_phc.update_file(hdphc_target_file, data)
    report_status(do_succeeded, "BH2")

    # -- HTOF -----
    limits = [-np.inf, np.inf]
    if f"{args.run_num:0=5}_{args.suffix}_htof" in phc_conf.limits_dict.keys():
        limits = phc_conf.limits_dict[f"{args.run_num:0=5}_{args.suffix}_htof"]
    data = update_phc.make_dictdata(os.path.join(script_dir, f"../results/root/run{args.run_num:0=5}_HTOF_PHC_{args.suffix}.root"), limits)
    do_succeeded = update_phc.update_file(hdphc_target_file, data)
    report_status(do_succeeded, "HTOF")


# elif args.suffix == "hdprm":
#     hdphc_dir = "{}/HDPHC".format(conf.param_dir)
#     hdphc_target_file = "{}/HodoPHC_run{:0=5}_{}".format(hdphc_dir, args.run_num, args.suffix)
#     if not os.path.isfile(hdphc_target_file):
#         shutil.copy("{}/HodoPHC_0".format(hdphc_dir), hdphc_target_file)
# ---------------------------------------------------------------------------




# # +--------+
# # | T0 ADC |
# # +--------+
# # -- up -----
# result1 = np.genfromtxt("{}/run{:0=5}_t0_ped_u_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# result2 = np.genfromtxt("{}/run{:0=5}_t0_mip_u_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# data = dict()
# for row1, row2 in zip(result1, result2):
#     # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
#     key = "2-0-{:.0f}-0-0".format(row1[0])
#     data[key] = [ row1[2], row2[2] ]
# update_file(data)
# print("update: T0 ADC up")



# # -- write USER and conf, runlist file  -----------------------------------
# def update_file(data, is_phc = False):
#     target_file = hdphc_target_file if is_phc else hdprm_target_file
#     key_length = 4 if is_phc else 5
#     buf = []
#     with open(target_file) as f:
#         for line in f:
#             s_list = line.split()
#             # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
#             # or
#             # CId - PlId - SegId - UorD(0:u, 1:d) - Type - nParam
#             key = s_list[0]
#             for i in range(1, key_length):
#                 key += "-"+s_list[i]
#             if key in data.keys():
#                 for i in range(len(data[key])):
#                     s_list[i+key_length] = data[key][i]                    
#             buf.append(s_list)

#     with open(target_file, mode='w') as f:
#         for l in buf:
#             f.write('\t'.join(str(item) for item in l))
#             f.write("\n")
# # ---------------------------------------------------------------------------

# print("\nbht_min: {}\nbht_max: {}\n".format(args.bht_min, args.bht_max)+"-"*20)
# # +---------+
# # | BHT TOT |
# # +---------+
# # -- up -----
# result = np.genfromtxt("{}/run{:0=5}_bht_tot_u_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# data = dict()
# for row in result:
#     # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
#     key = "1-0-{:.0f}-0-0".format(row[0])
#     data[key] = [ 0.0, row[2] ]
# update_file(data)
# print("update: BHT TOT up")

# # -- down -----
# result = np.genfromtxt("{}/run{:0=5}_bht_tot_d_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# data = dict()
# for row in result:
#     # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
#     key = "1-0-{:.0f}-0-1".format(row[0])
#     data[key] = [ 0.0, row[2] ]
# update_file(data)
# print("update: BHT TOT down")

# # +---------+
# # | BHT TDC |
# # +---------+
# # -- up -----
# result = np.genfromtxt("{}/run{:0=5}_bht_tdc_u_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# data = dict()
# tdc_mean = statistics.mean(result[ (args.bht_min <= result[:, 0])*(result[:, 0] <= args.bht_max) ][:, 2])
# for row in result:
#     # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
#     key = "1-0-{:.0f}-1-0".format(row[0])
#     if args.bht_min <= row[0] and row[0] <= args.bht_max:
#         data[key] = [ row[2] ]
#     else:
#         data[key] = [ np.round(tdc_mean, 1) ]
# update_file(data)
# print("update: BHT TDC up")

# # -- down -----
# result = np.genfromtxt("{}/run{:0=5}_bht_tdc_d_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# data = dict()
# tdc_mean = statistics.mean(result[ (args.bht_min <= result[:, 0])*(result[:, 0] <= args.bht_max) ][:, 2])
# for row in result:
#     # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
#     key = "1-0-{:.0f}-1-1".format(row[0])
#     if args.bht_min <= row[0] and row[0] <= args.bht_max:
#         data[key] = [ row[2] ]
#     else:
#         data[key] = [ np.round(tdc_mean, 1) ]
# update_file(data)
# print("update: BHT TDC down")


# # +--------+
# # | T0 ADC |
# # +--------+
# # -- up -----
# result1 = np.genfromtxt("{}/run{:0=5}_t0_ped_u_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# result2 = np.genfromtxt("{}/run{:0=5}_t0_mip_u_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# data = dict()
# for row1, row2 in zip(result1, result2):
#     # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
#     key = "2-0-{:.0f}-0-0".format(row1[0])
#     data[key] = [ row1[2], row2[2] ]
# update_file(data)
# print("update: T0 ADC up")

# # -- down -----
# result1 = np.genfromtxt("{}/run{:0=5}_t0_ped_d_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# result2 = np.genfromtxt("{}/run{:0=5}_t0_mip_d_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# data = dict()
# for row1, row2 in zip(result1, result2):
#     # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
#     key = "2-0-{:.0f}-0-1".format(row1[0])
#     data[key] = [ row1[2], row2[2] ]
# update_file(data)
# print("update: T0 ADC down")


# # +--------+
# # | T0 TDC |
# # +--------+
# # -- up -----
# result = np.genfromtxt("{}/run{:0=5}_t0_tdc_u_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# data = dict()
# for row in result:
#     key = "2-0-{:.0f}-1-0".format(row[0])
#     data[key] = [ row[2] ]
# update_file(data)
# print("update: T0 TDC up")

# # -- down -----
# result = np.genfromtxt("{}/run{:0=5}_t0_tdc_d_{}.csv".format(args.result_dir, args.run_num, args.suffix), delimiter=",", skip_header=1)
# data = dict()
# for row in result:
#     key = "2-0-{:.0f}-1-1".format(row[0])
#     data[key] = [ row[2] ]
# update_file(data)
# print("update: T0 TDC down")


# # +----------------+
# # | T0 Time Offset |
# # +----------------+
# t0_offset_path = "{}/run{:0=5}_t0_offset_{}.csv".format(args.result_dir, args.run_num, args.suffix)
# if os.path.isfile(t0_offset_path):
#     result = np.genfromtxt(t0_offset_path, delimiter=",", skip_header=1)
#     data = dict()
#     for row in result:
#         key = "2-0-{:.0f}-1-2".format(row[0])
#         data[key] = [ row[2], 1.0 ]
#     update_file(data)
#     print("update: T0 Time Offset")


# # +-----------------------------+
# # | BHT Pulse Height Correction |
# # +-----------------------------+
# # -- up -----
# bht_phc_path = "{}/run{:0=5}_bht_phc_u_{}.csv".format(args.result_dir, args.run_num, args.suffix)
# if os.path.isfile(bht_phc_path):
#     result = np.genfromtxt(bht_phc_path, delimiter=",", skip_header=1)
#     data = dict()
#     for row in result:
#         if args.bht_min <= row[0] and row[0] <= args.bht_max:
#             key = "1-0-{:.0f}-0".format(row[0])
#             data[key] = [ 1, 3, row[1], row[2], row[3] ]
#         else:
#             key = "1-0-{:.0f}-0".format(row[0])
#             data[key] = [ 0, 3, 0.0, 0.0, 0.0 ]
#     update_file(data, True)
#     print("update: BHT PHC up")


# # -- down -----
# bht_phc_path = "{}/run{:0=5}_bht_phc_d_{}.csv".format(args.result_dir, args.run_num, args.suffix)
# if os.path.isfile(bht_phc_path):
#     result = np.genfromtxt(bht_phc_path, delimiter=",", skip_header=1)
#     data = dict()
#     for row in result:
#         if args.bht_min <= row[0] and row[0] <= args.bht_max:
#             key = "1-0-{:.0f}-1".format(row[0])
#             data[key] = [ 1, 3, row[1], row[2], row[3] ]
#         else:
#             key = "1-0-{:.0f}-1".format(row[0])
#             data[key] = [ 0, 3, 0.0, 0.0, 0.0 ]
#     update_file(data, True)
#     print("update: BHT PHC down")


# # +----------------------------+
# # | T0 Pulse Height Correction |
# # +----------------------------+
# # -- up -----
# t0_phc_path = "{}/run{:0=5}_t0_phc_u_{}.csv".format(args.result_dir, args.run_num, args.suffix)
# if os.path.isfile(t0_phc_path):
#     result = np.genfromtxt(t0_phc_path, delimiter=",", skip_header=1)
#     data = dict()
#     for row in result:
#         key = "2-0-{:.0f}-0".format(row[0])
#         data[key] = [ 1, 3, row[1], row[2], row[3] ]
#     update_file(data, True)
#     print("update: T0 PHC up")

# # -- down -----
# t0_phc_path = "{}/run{:0=5}_t0_phc_d_{}.csv".format(args.result_dir, args.run_num, args.suffix)
# if os.path.isfile(t0_phc_path):
#     result = np.genfromtxt(t0_phc_path, delimiter=",", skip_header=1)
#     data = dict()
#     for row in result:
#         key = "2-0-{:.0f}-1".format(row[0])
#         data[key] = [ 1, 3, row[1], row[2], row[3] ]
#     update_file(data, True)
#     print("update: T0 PHC down")


# print("-"*20)