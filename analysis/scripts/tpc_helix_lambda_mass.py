#!/usr/bin/env python3
"""DstTPCHelixTracking 出力 ROOT（TTree ``tpc``、ブランチ ``lambda_mass``）から Λ 質量分布を描く。

依存: pip install uproot awkward matplotlib
使い方: python3 tpc_helix_lambda_mass.py path/to/file.root
"""

import sys

import awkward as ak
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import uproot
import os

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 28
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['axes.grid'] = False
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
plt.rcParams["xtick.major.size"] = 10                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 10                #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ

if len(sys.argv) != 2:
    sys.exit(f"usage: {sys.argv[0]} path/to/file.root")

ROOTFILE = sys.argv[1]

# --- 必要なら編集 ---
TREE = "tpc"
BINS = 400
XMIN, XMAX = 1.05, 1.25  # GeV/c^2（HistTools::BuildTPCHelixLambda と同じ）
# ------------------

with uproot.open(ROOTFILE) as src:
    if TREE not in src:
        sys.exit(f"error: tree {TREE!r} not in {ROOTFILE!r}")
    tree = src[TREE]
    if "lambda_mass" not in tree:
        sys.exit("error: branch 'lambda_mass' missing (EnableReconstructLambda ビルドの ROOT を使う)")
    jagged = tree["lambda_mass"].array(library="ak")

flat = ak.flatten(jagged, axis=None)
if len(flat) == 0:
    sys.exit("error: no lambda_mass entries")

masses = np.asarray(ak.to_numpy(flat), dtype=np.float64)

fig, ax = plt.subplots(figsize=(8, 7))
ax.hist(
    masses,
    bins=BINS,
    range=(XMIN, XMAX),
    color="k",
    histtype="step",
    zorder=2,
)
ax.text(
    0.5,
    0.5,
    "preliminary",
    transform=ax.transAxes,
    rotation=32,
    ha="center",
    va="center",
    fontsize=80,
    alpha=0.14,
    color="0.35",
    zorder=5,
)
ax.set_xlabel(r"$M(p\pi^-)$ [GeV/$c^2$]")
ax.set_ylabel("Counts")
ax.yaxis.set_major_formatter(mticker.EngFormatter(sep=""))
ax.set_xlim(1.07, 1.23)

plt.subplots_adjust(top = 0.98, bottom = 0.15, left = 0.15, right = 0.98, wspace=0, hspace=0)
script_dir = os.path.dirname(os.path.abspath(__file__))
img_path = os.path.join(script_dir, f"../../results/img/lambda_mass.png")
plt.savefig(img_path, format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()
