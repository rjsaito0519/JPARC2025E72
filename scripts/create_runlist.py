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
prefix = "dc" if args.dc else "hodo"
set1 = set(args.suffix)
set2 = set(["0", "Pi", "K"])

common_elements = set1.difference(set2)
if common_elements:
    print("please input correct suffix: 0, Pi, K")
    sys.exit()

# -- write conf file  -----------------------------------
for suffix in args.suffix:
    conf_dir = os.path.join(conf.analyzer_dir, "param", "conf")
    conf_target_file = f"{conf_dir}/analyzer_e72_run{args.run_num:0=5}_{prefix}_{suffix}.conf"
    if not os.path.isfile(conf_target_file):
        shutil.copy(os.path.join(conf_dir, "analyzer_e72_2025.conf"), conf_target_file)

    buf = []
    with open(conf_target_file) as f:
        for line in f:
            s_list = line.split()
            if len(s_list) != 0 and s_list[0] == "USER:":
                target = f"param/USER/UserParam_run{args.run_num:0=5}"
                if os.path.isfile(os.path.join(analyzer_dir, target)):
                    s_list[1] = target
            if prefix == "hodo":
                if len(s_list) != 0 and s_list[0] == "HDPRM:" and suffix in ["Pi", "K"]:
                    target = f"param/HDPRM/HodoParam_run{args.run_num:0=5}_{suffix}"
                    if os.path.isfile(os.path.join(analyzer_dir, target)):  
                        s_list[1] = target
                if len(s_list) != 0 and s_list[0] == "HDPHC:" and suffix in ["Pi", "K"]:
                    target = f"param/HDPHC/HodoPHC_run{args.run_num:0=5}_{suffix}"
                    if os.path.isfile(os.path.join(analyzer_dir, target)):  
                        s_list[1] = target
            # elif prefix == "dc":
            #     if len(s_list) != 0 and s_list[0] == "DCTDC:" and suffix in ["Pi_1", "K_1", "Pi_2", "K_2", "Pi_3", "K_3"]:
            #         s_list[1] = "param/DCTDC/DCTdcParam_run{:0=5}_{}".format(args.run_num, suffix[:suffix.find("_")])
            buf.append(s_list)

    with open(conf_target_file, mode='w') as f:
        for l in buf:
            f.write('\t'.join(str(item) for item in l))
            f.write("\n")
# ---------------------------------------------------------------------------

# -- write runlist file  -----------------------------------
runlist_dir = os.path.join(conf.analyzer_dir, "runmanager", "runlist")
runlist_target_file = os.path.join(runlist_dir, f"{prefix}_run{args.run_num:05d}.yml")
if not os.path.isfile(runlist_target_file):
    shutil.copy(os.path.join(runlist_dir, "myexample.yml"), runlist_target_file)

buf = []
with open(runlist_target_file) as f:
    for line in f:
        s_list = line.split()
        if len(s_list) != 0 and s_list[0] == "RUN:":
            buf.append(s_list)
            for suffix in args.suffix:
                buf.append([f"run{args.run_num:0=5}_{suffix}:"])                    
                if prefix == "hodo":
                    buf.append(["  bin: ./bin/Hodoscope"])
                elif prefix == "dc":
                    buf.append(["  bin: ./bin/BcOutTracking"])
                buf.append([f"  conf: ./param/conf/analyzer_e72_run{args.run_num:0=5}_{prefix}_{suffix}.conf"])
                buf.append([f"  data: ./rawdata/run{args.run_num:0=5}.dat"])
                buf.append([f"  root: ./root/{prefix}/{prefix}_run{args.run_num:0=5}_{suffix}.root"])
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
print("-"*30)
# ---------------------------------------------------------------------------