import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import numpy as np
import uproot
import os
import sys

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 28
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['axes.grid'] = True
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
plt.rcParams["xtick.major.size"] = 10                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 10                #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ
script_dir = os.path.dirname(os.path.abspath(__file__))


root_file_path = os.path.join(script_dir, "../data/root/run02293_COBO_diff.root")
file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

for it in tree["diff"]:
    print(np.mean(it)*9.765625e-04)


# fig = plt.figure(figsize=(8, 8))
# ax  = fig.add_subplot(111)



# x = tree[f"A{i}_spline_x"][0]
# y = tree[f"A{i}_spline_y"][0]        
# ax.plot(x, y)

# point_x = tree[f"A{i}_coeff_x"][0]
# point_y = tree[f"A{i}_coeff_y"][0]
# ax.plot(point_x, point_y, "o", ms = 7, color = "red")

# if (ax not in [ax3, ax6]):
#     ax.set_xticklabels([])
# ax.set_ylabel("$A_{}$".format(i), fontsize = 35)
# ax.yaxis.set_label_coords(-0.1, 0.5)

# fig.supxlabel(r"$p_{K}^{\rm Lab}$ [MeV/c]", x = 0.1 + (0.98-0.1)/2)
# fig.suptitle(arg_dict["title"], x = 0.1 + (0.98-0.1)/2)

# plt.subplots_adjust(left = 0.1, right=0.98, top=0.92, bottom = 0.14, hspace=0.02)
# img_save_path = os.path.join(script_dir, "../results/img/yield/example_{}.pdf".format(os.path.splitext( arg_dict["data"] )[0]))
# os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
# plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)

# plt.show()