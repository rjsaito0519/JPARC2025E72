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

def get_hist_data(file, key):
    hist_data = file[key].to_numpy()
    bin_values = hist_data[0]
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, bin_edges, bin_values

root_file_path = os.path.join(script_dir, "../data/root/run02293_COBO_diff.root")
file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

# for it in tree["diff"]:
#     print(np.mean(it)*9.765625e-04)


fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

center, edge, value = get_hist_data(file, "COBO_tdc_diff_2293_1")
ax1.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="C0", zorder = 3)
ax1.set_yscale("log")
ax1.xaxis.set_major_formatter(ptick.EngFormatter())
ax1.set_xlabel("TDC diff [ch]")

center, edge, value = get_hist_data(file, "COBO_time_diff_2293_1")
ax2.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="C0", zorder = 3)
ax2.set_yscale("log")
ax2.set_yticklabels([])
ax2.set_xlabel("Time diff [ns]")

plt.subplots_adjust(left = 0.1, right=0.98, top=0.92, bottom = 0.14, )
img_save_path = os.path.join(script_dir, "../results/img/COBO01.png")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)

plt.show()


fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)

center, edge, value = get_hist_data(file, "COBO_1st_hit_2293")
ax.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="C0", zorder = 3)
ax.set_yscale("log")
ax.xaxis.set_major_formatter(ptick.EngFormatter())
ax.set_xlabel("TDC [ch]")
ax.set_xlim(1.51E6, 1.64E6)

value = np.array(value)
edge  = np.array(edge)

threshold = 0.5 * np.max(value)

mask = value > threshold

x_min = edge[:-1][mask].min()
x_max = edge[1:][mask].max()

width = x_max - x_min
print(x_min, x_max)
print(f"Width (50% level): {width:.3e} ch")


plt.subplots_adjust(left = 0.1, right=0.98, top=0.92, bottom = 0.14, )
img_save_path = os.path.join(script_dir, "../results/img/COBO02.png")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)

plt.show()