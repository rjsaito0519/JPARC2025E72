#!/usr/bin/env python3

# ---------------------------------------------------------------------------
import argparse
from curses.ascii import SUB
parser = argparse.ArgumentParser(
    prog="set_tdc_window",
    usage="python3 set_tdc_window.py <rootfile_path> <counter_name>",
    description="set tdc window",
    epilog="end",
    add_help=True,
)
parser.add_argument("rootfile_path", type=str, help="Input rootfile path")
parser.add_argument("counter_name", type=str, help="Input counter name")
parser.add_argument('--check', action="store_true", help='Set check to True')
args = parser.parse_args()
# ---------------------------------------------------------------------------

import uproot
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

# ---------------------------------------------------------------------------
import conf
filename = os.path.basename(args.rootfile_path)
index = filename.find("run")
run_num = int(filename[index+3:index+8])

SUB_DIR = "e72"

script_dir = os.path.dirname(os.path.abspath(__file__))
# user_param_dir = os.path.join(conf.analyzer_dir, "param", "USER", SUB_DIR)
user_param_dir = os.path.join(script_dir, "param", "USER", SUB_DIR)

target_file = os.path.join(user_param_dir, f"UserParam_run{run_num:05d}")

if not os.path.isfile(target_file):
    shutil.copy(os.path.join(user_param_dir, "UserParam_e72_20250518"), target_file)
# ---------------------------------------------------------------------------

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 22
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
plt.rcParams['figure.subplot.left'] = 0.1
plt.rcParams['figure.subplot.right'] = 0.98
plt.rcParams['figure.subplot.top'] = 0.95
plt.rcParams['figure.subplot.bottom'] = 0.1

num_of_ch = dict(
    BHT = 63,
    T0 = 5,
    BAC = 4,
    SAC3 = 1,
    KVC = 8,
    BH2 = 1,
    HTOF = 34,
    T1 = 1,
    CVC = 8,
    SFV = 1
)

# matplotlibのグラフを動的にする
# ---------------------------------------------------------------------------
fill_patch = None
vline_patch = None
bht_mode = False
selected_points = []
def on_click(event):
    global fill_patch
    global vline_patch
    global bht_mode

    # 右クリック（button=3）かどうかをチェック
    if event.button == 3:
        # クリックされた点のx, yデータを取得
        x = event.xdata
        print(f'Right click at x={x}')
        
        # 選択した点をリストに追加
        selected_points.append(x)
        
        # 選択した点が2つになったら塗りつぶす
        if len(selected_points) > 1:
            if fill_patch:
                fill_patch.remove()
            if vline_patch:
                vline_patch.remove()
            x1 = selected_points[-2]
            x2 = selected_points[-1]
            x3 = np.nan
            if bht_mode:
                x1 = selected_points[-3]
                x2 = selected_points[-2]
                x3 = selected_points[-1]

            # xの範囲を決定
            x_min, x_max = sorted([x1, x2])
            # yの範囲を決定
            y_min, y_max = ax.get_ylim()
            
            # 塗りつぶし
            fill_patch = ax.fill_betweenx([y_min, y_max], x_min, x_max, color='gray', alpha=0.5)
            if bht_mode:
                vline_patch = ax.vlines(x3, y_min, y_max, color='red', ls = "dashed")
            plt.draw()

# キーボードイベントのハンドラ関数
def on_key(event):
    global bht_mode          
    if event.key == "w":
        if len(selected_points) > 1:

            tdc_label = "{}_TDC".format(args.counter_name)

            buf = []
            do_update = False
            log = []
            with open(target_file) as f:
                for line in f:
                    s_list = line.split()                
                    if len(s_list) != 0:
                        if bht_mode:
                            x1 = selected_points[-3]
                            x2 = selected_points[-2]
                            x3 = selected_points[-1]
                            x_min, x_max = sorted([x1, x2])
                            if s_list[0] == tdc_label:
                                s_list[1] = int(x_min+0.5)
                                s_list[2] = int(x_max+0.5)  
                                log.append(s_list)                              
                            elif s_list[0] == tdc_label+"_Pi":
                                s_list[1] = int(x_min+0.5)
                                s_list[2] = int(x3+0.5)
                                log.append(s_list)                   
                            elif s_list[0] == tdc_label+"_K":
                                s_list[1] = int(x3+0.5)
                                s_list[2] = int(x_max+0.5)  
                                log.append(s_list)
                                do_update = True
                        else:
                            if s_list[0] == tdc_label:
                                x1 = selected_points[-2]
                                x2 = selected_points[-1]
                                x_min, x_max = sorted([x1, x2])
                                s_list[1] = int(x_min+0.5)
                                s_list[2] = int(x_max+0.5)
                                log.append(s_list)
                                do_update = True
                    buf.append(s_list)

            with open(target_file, mode='w') as f:
                for l in buf:
                    f.write('\t'.join(str(item) for item in l))
                    f.write("\n")

            if do_update:
                print("==========\nupdate " + tdc_label)
                for line in log:
                    print(line)
                print("==========")
                
    elif event.key == "1":
        print("bht mode ON")
        bht_mode = True
    elif event.key == "2":
        print("bht mode OFF")
        bht_mode = False
# ---------------------------------------------------------------------------

file = uproot.open(args.rootfile_path)

if args.counter_name in ["BHT", "T0", "T1", "CVC", "NC", "BH2", "HTOF", "SAC3", "SFV"]:
    hist_data = file["{}_TDC_seg0U".format(args.counter_name)].to_numpy()

    bin_values_u = np.zeros_like(hist_data[0])
    bin_values_d = np.zeros_like(hist_data[0])
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    for i in range(num_of_ch[args.counter_name]):
        hist_data_u = file["{}_TDC_seg{}U".format(args.counter_name, i)].to_numpy()
        hist_data_d = file["{}_TDC_seg{}D".format(args.counter_name, i)].to_numpy()
        bin_values_u += hist_data_u[0]
        bin_values_d += hist_data_d[0]
        
    fig = plt.figure(figsize=(8, 8))
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)
    ax = fig.add_subplot(111)

    # プロット
    ax.hist(bin_centers, bins=bin_edges, weights=bin_values_u, histtype='step', color="C1", label = "up")
    ax.hist(bin_centers, bins=bin_edges, weights=bin_values_d, histtype='step', color="C0", label = "down")
    ax.set_yscale("log")
    
    if (args.check):
        tdc_label = "{}_TDC".format(args.counter_name)
        with open(target_file) as f:
            for line in f:
                s_list = line.split()                
                if len(s_list) != 0:
                    if s_list[0] == tdc_label:
                        x_min = float(s_list[1])
                        x_max = float(s_list[2])
                    elif s_list[0] == tdc_label+"_Pi":
                        x_sep = float(s_list[2])

        y_min = 0
        y_max = np.max([np.max(bin_values_u), np.max(bin_values_d)])
        ax.fill_betweenx([y_min, y_max], x_min, x_max, color='gray', alpha=0.25)
        if args.counter_name == "BHT":
            ax.vlines(x_sep, y_min, y_max, color='red', ls = "dashed")
        plt.legend()

    plt.show()

elif args.counter_name in ["BAC", "SAC"]:
    hist_data = file["{}_TDC_seg{}U".format(args.counter_name, num_of_ch[args.counter_name])].to_numpy()

    bin_values = hist_data[0]
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    fig = plt.figure(figsize=(8, 8))
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)
    ax = fig.add_subplot(111)

    # プロット
    ax.hist(bin_centers, bins=bin_edges, weights=bin_values, histtype='step')
    ax.set_yscale("log")
    
    if (args.check):
        tdc_label = "{}_TDC".format(args.counter_name)
        with open(target_file) as f:
            for line in f:
                s_list = line.split()                
                if len(s_list) != 0:
                    if s_list[0] == tdc_label:
                        x_min = float(s_list[1])
                        x_max = float(s_list[2])

        y_min = 0
        y_max = np.max(bin_values)
        ax.fill_betweenx([y_min, y_max], x_min, x_max, color='gray', alpha=0.25)

    plt.show()

elif args.counter_name in ["KVC", "KVC1"]:
    hist_data = file["{}_TDC_seg0S".format(args.counter_name)].to_numpy()

    bin_values = np.zeros_like(hist_data[0])
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    for i in range(num_of_ch[args.counter_name]):
        hist_data = file["{}_TDC_seg{}S".format(args.counter_name, i)].to_numpy()
        bin_values += hist_data[0]
        
    fig = plt.figure(figsize=(8, 8))
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)
    ax = fig.add_subplot(111)

    # プロット
    ax.hist(bin_centers, bins=bin_edges, weights=bin_values, histtype='step', color="C0")
    ax.set_yscale("log")
    
    if (args.check):
        tdc_label = "{}_TDC".format(args.counter_name)
        with open(target_file) as f:
            for line in f:
                s_list = line.split()                
                if len(s_list) != 0:
                    if s_list[0] == tdc_label:
                        x_min = float(s_list[1])
                        x_max = float(s_list[2])

        y_min = 0
        y_max = np.max(bin_values)
        ax.fill_betweenx([y_min, y_max], x_min, x_max, color='gray', alpha=0.25)

    plt.show()
