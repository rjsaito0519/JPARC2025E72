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
