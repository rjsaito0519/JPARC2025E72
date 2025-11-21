#!/usr/bin/env python3

# ---------------------------------------------------------------------------
import argparse
parser = argparse.ArgumentParser(
    prog="make_runlist",
    usage="python3 create_runlist.py <run_num> <suffix> <--dc>",
    description="",
    epilog="end", 
    add_help=True,
)
parser.add_argument("run_num", type=int, help="Input run number")
parser.add_argument("suffix",  type=str, help="Input suffix", nargs='*')
parser.add_argument('--dc', action="store_true", help='Set check to True')
args = parser.parse_args()
# ---------------------------------------------------------------------------

import sys
import os
import shutil
import conf
SUB_DIR = "e72"
for param_type in ["conf", "USER", "HDPRM", "HDPHC", "DCTDC", "DCDRFT"]:
    target_dir = os.path.join(conf.analyzer_dir, "param", SUB_DIR, param_type)
    os.makedirs(target_dir, exist_ok=True)

prefix = "dc" if args.dc else "hodo"
set1 = set(args.suffix)
set2 = set(["0", "Pi_hdprm", "K_hdprm", "Pi_t0", "K_t0", "Pi_hdphc", "K_hdphc"])
set3 = set(["0", "Pi_tdc", "K_tdc", "Pi_drift", "K_drift"])

if prefix == "hodo":
    common_elements = set1.difference(set2)
    if common_elements or len(set1) == 0:
        print("please input correct suffix: 0, Pi/K_hdprm, Pi/K_t0, Pi/K_hdphc")
        sys.exit()
elif prefix == "dc":
    common_elements = set1.difference(set3)
    if common_elements or len(set1) == 0:
        print("please input correct suffix: 0, Pi/K_tdc, Pi/K_drift")
        sys.exit()

# -- write conf file  -----------------------------------
for suffix in args.suffix:
    suffix_head = suffix.split("_")[0]
    conf_dir = os.path.join(conf.analyzer_dir, "param", "conf")
    conf_target_file = f"{conf_dir}/{SUB_DIR}/analyzer_run{args.run_num:0=5}_{prefix}_{suffix_head}.conf"
    if not os.path.isfile(conf_target_file):
        shutil.copy(os.path.join(conf_dir, "analyzer_e72_20251113.conf"), conf_target_file)

    buf = []
    with open(conf_target_file) as f:
        for line in f:
            s_list = line.split()
            if len(s_list) != 0 and s_list[0] == "USER:":
                target = f"param/USER/{SUB_DIR}/UserParam_run{args.run_num:0=5}"
                if not os.path.isfile(os.path.join(conf.analyzer_dir, target)):
                    shutil.copy(os.path.join(conf.analyzer_dir, "param/USER/UserParam_e72_20251104"), os.path.join(conf.analyzer_dir, target))
                s_list[1] = target
            if len(s_list) != 0 and s_list[0] == "HDPRM:":
                target = f"param/HDPRM/{SUB_DIR}/HodoParam_run{args.run_num:0=5}_{suffix_head}"
                if os.path.isfile(os.path.join(conf.analyzer_dir, target)):  
                    s_list[1] = target
            if len(s_list) != 0 and s_list[0] == "HDPHC:":
                target = f"param/HDPHC/{SUB_DIR}/HodoPHCParam_run{args.run_num:0=5}_{suffix_head}"
                if os.path.isfile(os.path.join(conf.analyzer_dir, target)):  
                    s_list[1] = target
            if len(s_list) != 0 and s_list[0] == "DCTDC:":
                target = f"param/DCTDC/{SUB_DIR}/DCTdcParam_run{args.run_num:0=5}_{suffix_head}"
                if os.path.isfile(os.path.join(conf.analyzer_dir, target)):  
                    s_list[1] = target
            buf.append(s_list)

    with open(conf_target_file, mode='w') as f:
        for l in buf:
            f.write('\t'.join(str(item) for item in l))
            f.write("\n")
# ---------------------------------------------------------------------------

# -- write runlist file  -----------------------------------
runlist_dir = os.path.join(conf.analyzer_dir, "runmanager", "runlist")
runlist_target_file = os.path.join(runlist_dir, f"{SUB_DIR}/{prefix}_run{args.run_num:05d}.yml")
if not os.path.isfile(runlist_target_file):
    shutil.copy(os.path.join(runlist_dir, "myexample.yml"), runlist_target_file)

buf = []
with open(runlist_target_file) as f:
    for line in f:
        s_list = line.split()
        if len(s_list) != 0 and s_list[0] == "RUN:":
            buf.append(s_list)
            for suffix in args.suffix:
                buf.append([f"run{args.run_num:0=5}_{suffix}_in:"])
                if prefix == "hodo":
                    buf.append(["  bin: ./bin/Hodoscope"])
                    buf.append([f"  conf: ./param/conf/{SUB_DIR}/analyzer_run{args.run_num:0=5}_{prefix}_{suffix_head}.conf"])
                    buf.append([f"  data: ./rawdata/run{args.run_num:0=5}.dat"])
                    buf.append(["  root: {}".format(os.path.join(conf.output_dir, f"{prefix}_run{args.run_num:0=5}_{suffix}.root"))])
                elif prefix == "dc":
                    buf.append(["  bin: ./bin/BcInTracking"])
                    buf.append([f"  conf: ./param/conf/{SUB_DIR}/analyzer_run{args.run_num:0=5}_{prefix}_{suffix_head}.conf"])
                    buf.append([f"  data: ./rawdata/run{args.run_num:0=5}.dat"])
                    buf.append(["  root: {}".format(os.path.join(conf.output_dir, f"{prefix}_in_run{args.run_num:0=5}_{suffix}.root"))])
                    buf.append([f"run{args.run_num:0=5}_{suffix}_out:"])
                    buf.append(["  bin: ./bin/BcOutTracking"])
                    buf.append([f"  conf: ./param/conf/{SUB_DIR}/analyzer_run{args.run_num:0=5}_{prefix}_{suffix_head}.conf"])
                    buf.append([f"  data: ./rawdata/run{args.run_num:0=5}.dat"])
                    buf.append(["  root: {}".format(os.path.join(conf.output_dir, f"{prefix}_out_run{args.run_num:0=5}_{suffix}.root"))])
            break
        else:
            buf.append(s_list)

with open(runlist_target_file, mode='w') as f:
    for l in buf:
        if len(l) > 0 and (l[0][0] == "#" or l[0] in ["WORKDIR:", "DEFAULT:", "RUN:"]):
            f.write(' '.join(str(item) for item in l))
            f.write("\n")
        else:
            f.write("  ")
            f.write(' '.join(str(item) for item in l))
            f.write("\n")
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
print("\n" + "-"*30)
print("make runlist: {}_run{:0=5}.yml".format(prefix, args.run_num))
with open(runlist_target_file, "r") as f:
    print(f.read())
print("-"*30)
# ---------------------------------------------------------------------------