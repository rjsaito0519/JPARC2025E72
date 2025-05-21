import uproot
import os 

script_dir = os.path.dirname(os.path.abspath(__file__))
root_file_path = os.path.join(script_dir, "../results/root/run00088_T0_HDPEM_Pi.root")

if "T0" in root_file_path:
    print("aaaaa")

file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

print(tree["tdc_p0_val"])

data = dict()
for i in range(len(tree["adc_p0_val"])):
    # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
    for UorD in [0, 1]:
        # -- ADC -----
        key = f"2-0-{i:.0f}-0-{UorD:.0f}"
        data[key] = [ tree["adc_p0_val"][i][UorD], tree["adc_p1_val"][i][UorD] ]

        # -- TDC -----
        key = f"2-0-{i:.0f}-1-{UorD:.0f}"
        data[key] = [ tree["tdc_p0_val"][i][UorD], -0.0009765625 ]

# print(data)

# for row1, row2 in zip(result1, result2):
#     # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
#     key = "2-0-{:.0f}-0-0".format(row1[0])
#     data[key] = [ row1[2], row2[2] ]
# update_file(data)
# print("update: T0 ADC up")

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