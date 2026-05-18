#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DstTPCHelixTracking 出力（ROOT ツリー ``tpc``）用の 3D イベントディスプレイ（matplotlib）。

可視化のメモ（ROOT / matplotlib 共通の並べ方）:
  **プロットの第1・第2・第3座標に渡す値の並び = 局所 (z, x, y)**。すなわち
  ``(Xv, Yv, Zv) = (z_{\\mathrm{TPC}}, x_{\\mathrm{TPC}}, y_{\\mathrm{TPC}})``。
  ツリーの ``hitpos_x`` 等は局所 ``x,y,z`` なので、表示では ``plot(z, x, y)`` 相当になる。

前提:
  - 螺旋パラメータ ``helix_cx``, …, ``helix_t`` と測定座標 ``hitpos_x/y/z`` を使用する。
  - （推奨）Dst が ``helix_theta_min`` / ``helix_theta_max``（または旧 ``helix_t_min`` / ``helix_t_max``）
    を書き出す場合は、フィッタの θ 範囲で螺旋線を描き、ヒット列からの推定より安定する。
  - ツリー上の螺旋パラメータは **フィッティング用内部座標**（``TPCLocalTrackHelix::LocalPosition``）。
    表示は ``LocalToGlobal`` 後の **TPC 局所座標**（``hitpos_*`` と同じ系）に合わせる::

      xi = cx + r*cos(t),  yi = cy + r*sin(t),  zi = z0 + dz*r*t
      x = -xi,  y = zi,  z = yi + Z_TARGET   （Z_TARGET = -143 mm, ``tpc::Z_TARGET``）

  - matplotlib の 3D 表示は track 用マクロと同じ可視化回転（上記の (z,x,y) 順）を掛ける。
    ターゲット円筒は **局所 y が母線**、断面は **局所 x-z**（中心 z≈-143 mm）。幾何は ``check_tpc_helix_track_3d.C`` / 本 Python で Helix 用に維持する（``check_tpc_track_3d.C`` 本体は変更しない）。

PID（dE/dx 候補のビットパターン; 物理種の確定ではない）:
  - 0x1: pion  → 表示色 緑
  - 0x2: kaon  → 表示色 青
  - 0x4: proton → 表示色 赤
  複数ビットが立つ場合は表示色の優先順位を **Kaon > Proton > Pion** に固定する
  （コンソールには候補文字列も併記する）。

注意:
  - 螺旋線の θ 区間は ``check_tpc_helix_track_3d.C`` の ``HelixThetaDrawRange`` と同期:
    2π 超の広がり・先頭末尾との不整合・6–94% タイルより異常に広い区間では先頭〜末尾または
    タイル幅に縮め、端にわずかなマージンを付ける。θ 区間に応じてサンプル数を増やす。
    有効な θ 区間が得られ、``helix_r`` が有限なら螺旋線を描く（半径 mm が小さくてもスキップしない）。
  - DstTPCHelixTracking の既定ビルドでは RawHit/RawCluster が無効なことが多く、
    生ヒット・全クラスタ枝は存在しない。螺旋＋トラック上ヒット中心の表示となる。
  - ビルドフラグにより ``calpos_*``, ``vtxTpc`` 等の枝が欠ける場合がある。その場合は
    該当描画をスキップする。
  - ``isBeam`` は Dst 側の実装状況によっては **PyROOT が読むとセグフォ** することがある。
    既定では **一切読み込まない**。安定して枝がある環境では
    ``set_path_helix_3d(path, read_isbeam=True)`` または ``--read-isbeam`` を付ける。
  - **既定 I/O は uproot**（PyROOT/cppyy を避ける）。``pip install uproot awkward`` が必要。
    どうしても PyROOT 使う場合のみ ``--pyroot`` または ``backend=\"pyroot\"``。

Usage:
  python check_tpc_helix_track_3d.py <rootfile> [entry]
  python check_tpc_helix_track_3d.py <rootfile> --run 12345
  python check_tpc_helix_track_3d.py <rootfile> --run 12345 --event 678
      # run_number（と任意で event_number）が一致する最初のエントリ
  python check_tpc_helix_track_3d.py <rootfile> --event 678
      # event_number のみ一致（先頭から最初の 1 件。複数ランで重複しうるので注意）
  python check_tpc_helix_track_3d.py <rootfile> --kp
      # dE/dx PID ビットで Kaon(0x2) または Proton(0x4) を含むトラックだけ描画
  python check_tpc_helix_track_3d.py <rootfile> [entry] --save
      # 透過背景・高 DPI の PNG（タイトル・凡例なし）。ファイル名省略時は helix_run*_ev*_entry*.png
  python check_tpc_helix_track_3d.py <rootfile> --event 123 --save fig.png --save-dpi 400
  python check_tpc_helix_track_3d.py <rootfile> --browse
      # 図ウィンドウをフォーカスしてキー「n」で次のツリーエントリ。閉じるまで終了
  python check_tpc_helix_track_3d.py <rootfile> --scan-kp --close-dist-max 20 [--max-scan N]
  (--scan-kp 単独は不可。isBeam を読まないときは --close-dist-max 必須)

  python
  >>> from check_tpc_helix_track_3d import *
  >>> set_path_helix_3d("path/to/file.root")  # uproot 既定（pip install uproot awkward）
  >>> set_path_helix_3d("path/to/file.root", read_isbeam=True)  # isBeam を読む場合
  >>> set_path_helix_3d("path/to/file.root", backend="pyroot")  # 旧来の PyROOT（非推奨）
  >>> event_helix_3d()              # random entry
  >>> event_helix_3d(123)         # specific entry
  >>> event_helix_3d(123, kp_only=True)  # K/p ビットのトラックだけ描画
  >>> event_helix_3d(123, save_png="out.png")  # 透過 PNG（タイトル・凡例なし）
  >>> event_helix_3d(123, save_png="", save_dpi=400)  # 自動ファイル名
  >>> event_helix_3d_browse()     # 図上で「n」で次エントリ（ウィンドウを閉じて終了）
  >>> event_helix_3d_browse(123)  # 最初は entry 123、以降は図上で「n」
  >>> e = find_entry_helix_3d(12345, 678)  # Run/Event に対応するツリーエントリを探索
  >>> if e is not None: event_helix_3d(e)
  >>> e = scan_kp_candidate(0, 5000, close_dist_max=15.0)
  >>> if e is not None: event_helix_3d(e)

Controls:
  - 図ウィンドウにフォーカスがあるとき キー ``n``: 次のツリーエントリ（末尾の次は先頭へループ）
  - Left mouse + drag: Rotate
  - Right mouse + drag: Zoom
  - Middle mouse + drag: Pan
  - Mouse wheel: Zoom in/out
  - 既定では matplotlib の 3D 軸枠は使わず、角の小さな矢印のみ（ラベルは局所 z,x,y＝表示の第1〜3軸方向）
"""

from __future__ import annotations

import argparse
import math
import os
import random
import sys
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np

# tpc::Z_TARGET (include/TPCPadHelper.hh) — LocalToGlobal の z オフセット
TPC_Z_TARGET = -143.0

# check_tpc_track_3d.C の SetRange と同じ表示空間 (局所 x,y,z)→(表示 X,Y,Z)=(z,x,y) [mm]
_XMIN_O, _XMAX_O = -311.0, 311.0
_YMIN_O, _YMAX_O = -310.0, 310.0
_ZMIN_O, _ZMAX_O = -311.0, 311.0
DISPLAY_XLIM = (_ZMIN_O, _ZMAX_O)
DISPLAY_YLIM = (_XMIN_O, _XMAX_O)
DISPLAY_ZLIM = (_YMIN_O, _YMAX_O)

# 3D 初期カメラ (Axes3D.view_init)。表示 (matplotlib X,Y,Z) = (局所 z, 局所 x, 局所 y)。
#
# 要件: 視線を局所 x に平行、カメラは局所 -x 側（＝ matplotlib の Y が負の側から +Y 向きに見る）。
# matplotlib では「+Y 軸方向に視線を向ける」＝ XZ 平面を正面にする横から見た形で、
# elev=0 がその平面内（局所 z–局所 y 平面）の視点。azim の符号で -x / +x のどちら側かが入れ替わる。
# 画角が反対なら HELIX_VIEW_AZIM の符号を反転（例: -90 <-> 90）。
# elev を少し上げると「上から覗く」感が出る（大きすぎると横視点から離れる）。
HELIX_VIEW_ELEV = 50.0
HELIX_VIEW_AZIM = -90.0

# --save 用 PNG（既定 DPI。fig は保存時のみやや大きめにする）
HELIX_SAVE_DPI = 300
HELIX_SAVE_FIGSIZE_IN = (12.0, 12.0)

# --- globals ---
g_backend = ""  # "uproot" | "pyroot"
g_file_helix = None
g_tree_helix = None
g_uproot_file = None
g_uproot_tree = None
g_entries_helix = 0
g_fig_helix = None
g_ax_helix = None
g_current_helix_entry: Optional[int] = None
g_helix_n_samples: int = 160
g_helix_kp_only: bool = False  # True のとき K(0x2) または p(0x4) ビットを持つトラックのみ描画
g_helix_png_export: bool = False  # True のときタイトル・凡例なし（--save）

# optional branch presence (filled in set_path_helix_3d)
g_has_calpos = False
g_has_vertex = False
g_has_close_dist = False
g_has_helix_t_range = False  # helix_theta_min/max または helix_t_min/max のいずれかがツリーにある
g_helix_trange_src = 0  # 0=なし, 1=helix_theta_min/max（現行 Dst）, 2=helix_t_min/max（旧名）
# isBeam: branch may exist but still crash if not wired; only read when read_isbeam=True
g_branch_isbeam = False
g_read_isbeam = False


def _close_helix_sources() -> None:
    global g_file_helix, g_tree_helix, g_uproot_file, g_uproot_tree, g_entries_helix, g_backend
    global g_has_helix_t_range, g_helix_trange_src
    if g_uproot_file is not None:
        try:
            g_uproot_file.close()
        except Exception:
            pass
        g_uproot_file = None
        g_uproot_tree = None
    if g_file_helix is not None:
        try:
            g_file_helix.Close()
        except Exception:
            pass
        g_file_helix = None
    g_tree_helix = None
    g_entries_helix = 0
    g_backend = ""
    g_has_helix_t_range = False
    g_helix_trange_src = 0


def _io_ready() -> bool:
    if g_entries_helix <= 0:
        return False
    if g_backend == "uproot":
        return g_uproot_tree is not None
    if g_backend == "pyroot":
        return g_tree_helix is not None
    return False


def _uproot_has_branch(tree, name: str) -> bool:
    try:
        return name in tree.keys()
    except Exception:
        return False


def _ak_to_py(x):
    """awkward の 1 イベント分の値を純粋な Python に。"""
    import awkward as ak

    if isinstance(x, (int, float, bool)):
        return x
    try:
        return ak.to_list(x)
    except Exception:
        return x


def _cell0_scalar_int(chunk, name: str) -> int:
    v = chunk[name][0]
    if hasattr(v, "item"):
        return int(v.item())
    return int(v)


def _cell0_vec_int(chunk, name: str) -> List[int]:
    raw = _ak_to_py(chunk[name][0])
    if raw is None:
        return []
    if isinstance(raw, (int, float)):
        return [int(raw)]
    return [int(z) for z in raw]


def _cell0_vec_float(chunk, name: str) -> List[float]:
    raw = _ak_to_py(chunk[name][0])
    if raw is None:
        return []
    if isinstance(raw, (int, float)):
        return [float(raw)]
    return [float(z) for z in raw]


def _cell0_nested_float(chunk, name: str, n_tracks: int) -> List[List[float]]:
    raw = _ak_to_py(chunk[name][0])
    if raw is None:
        return [[] for _ in range(max(n_tracks, 0))]
    if n_tracks <= 0:
        return []
    if not isinstance(raw, (list, tuple)):
        return [[float(raw)]]
    if not raw:
        return [[] for _ in range(n_tracks)]
    out: List[List[float]] = []
    for i in range(n_tracks):
        if i < len(raw) and isinstance(raw[i], (list, tuple)):
            out.append([float(z) for z in raw[i]])
        elif i < len(raw):
            out.append([float(raw[i])])
        else:
            out.append([])
    return out


def _ensure_pyroot():
    """PyROOT を遅延 import。matplotlib 等との干渉を少し抑える。"""
    import ROOT

    try:
        ROOT.PyConfig.IgnoreCommandlineOptions = True
    except Exception:
        pass
    try:
        ROOT.gROOT.SetBatch(True)
    except Exception:
        pass
    return ROOT


@dataclass
class EventDataHelix3D:
    runnum: int = 0
    evnum: int = 0
    ntTpc: int = 0
    nhtrack: List[int] = field(default_factory=list)
    is_beam: List[int] = field(default_factory=list)
    charge: List[int] = field(default_factory=list)
    mom0: List[float] = field(default_factory=list)
    pid: List[int] = field(default_factory=list)
    helix_cx: List[float] = field(default_factory=list)
    helix_cy: List[float] = field(default_factory=list)
    helix_z0: List[float] = field(default_factory=list)
    helix_r: List[float] = field(default_factory=list)
    helix_dz: List[float] = field(default_factory=list)
    helix_t_min: List[float] = field(default_factory=list)
    helix_t_max: List[float] = field(default_factory=list)
    helix_t: List[List[float]] = field(default_factory=list)
    hitpos_x: List[List[float]] = field(default_factory=list)
    hitpos_y: List[List[float]] = field(default_factory=list)
    hitpos_z: List[List[float]] = field(default_factory=list)
    calpos_x: List[List[float]] = field(default_factory=list)
    calpos_y: List[List[float]] = field(default_factory=list)
    calpos_z: List[List[float]] = field(default_factory=list)
    vtxTpc: List[List[float]] = field(default_factory=list)
    vtyTpc: List[List[float]] = field(default_factory=list)
    vtzTpc: List[List[float]] = field(default_factory=list)
    closeDistTpc: List[List[float]] = field(default_factory=list)


g_event_helix = EventDataHelix3D()


def _vec_size(obj) -> int:
    """PyROOT の std::vector 互換の長さ。"""
    try:
        return int(obj.size())
    except Exception:
        try:
            return len(obj)
        except Exception:
            return 0


def _vec_int_to_list(obj) -> List[int]:
    n = _vec_size(obj)
    return [int(obj[i]) for i in range(n)]


def _vec_double_to_list(obj) -> List[float]:
    n = _vec_size(obj)
    return [float(obj[i]) for i in range(n)]


def _nested_vec_double_from_tree(tree, attr: str, n_tracks: int) -> List[List[float]]:
    """std::vector<std::vector<double>> を要素単位でコピー（list(枝) より cppyy で安定しやすい）。"""
    obj = getattr(tree, attr)
    out: List[List[float]] = []
    for i in range(n_tracks):
        if i < _vec_size(obj):
            row = obj[i]
            out.append([float(row[j]) for j in range(_vec_size(row))])
        else:
            out.append([])
    return out


def decode_pid_candidates(pid_code: int) -> str:
    """TPCEventAnalyzer::DecodePidCandidates と同じ規則（pid==0 は電子候補扱い）。"""
    if pid_code == 0:
        return "e"
    names = ""
    if pid_code & 0x1:
        names += "pi"
    if pid_code & 0x2:
        if names:
            names += "+"
        names += "K"
    if pid_code & 0x4:
        if names:
            names += "+"
        names += "p"
    return names if names else "none"


def track_visible_kp_filter(pid_code: int) -> bool:
    """
    ``g_helix_kp_only`` が False なら常に True。
    True のときは PID ビットに Kaon(0x2) または Proton(0x4) が含まれるトラックだけ True。
    """
    if not g_helix_kp_only:
        return True
    return (pid_code & 0x2) != 0 or (pid_code & 0x4) != 0


def pid_display_color(pid_code: int) -> Tuple[str, str]:
    """
    (matplotlib color, label key).
    複数ビット: Kaon > Proton > Pion の優先。
    """
    if pid_code & 0x2:
        return "tab:blue", "K"
    if pid_code & 0x4:
        return "tab:red", "p"
    if pid_code & 0x1:
        return "tab:green", "pi"
    return "0.45", "?"  # gray-like


def helix_xyz_internal(cx: float, cy: float, z0: float, r: float, dz: float, theta: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    フィッタ内部座標 (xi, yi, zi) = ``TPCLocalTrackHelix.cc`` の ``LocalPosition(par,t)`` と同じ::

        xi = cx + r*cos(t),  yi = cy + r*sin(t),  zi = z0 + dz*r*t   （C++ は par[4]*par[3]*t）
    """
    xi = cx + r * np.cos(theta)
    yi = cy + r * np.sin(theta)
    zi = z0 + dz * r * theta
    return xi, yi, zi


def helix_xyz(cx: float, cy: float, z0: float, r: float, dz: float, theta: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    TPC 局所 (x,y,z) = ``LocalToGlobal(LocalPosition)``（``TPCLocalTrackHelix.cc``）=
    ``TPCLTrackHit::GetLocalCalPosHelix`` と同じ系::

        x = -xi,  y = zi,  z = yi + Z_TARGET

    ツリーの ``hitpos_*`` は ``GetLocalHitPos()``、螺旋上の参照点は ``calpos_*`` = ``GetLocalCalPosHelix()``。
    螺旋**線**は上式（cal と同じ）で描くので、hit と cal の差（残差）ぶれは残る。
    """
    xi, yi, zi = helix_xyz_internal(cx, cy, z0, r, dz, theta)
    x = -xi
    y = zi
    z = yi + TPC_Z_TARGET
    return x, y, z


def tpc_local_to_display(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """check_tpc_track_3d.C と同じ (局所 x,y,z) → (表示 X,Y,Z) = (z, x, y)。"""
    return float(z), float(x), float(y)


def tpc_local_to_display_vec(
    x: np.ndarray, y: np.ndarray, z: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    return z, x, y


def _make_tpc_frame_segment_arrays() -> Tuple[np.ndarray, np.ndarray]:
    """
    check_tpc_track_3d.C の八角柱＋ターゲット円筒をワイヤー線分 (N,2,3) にする。
    座標は表示空間 (z,x,y)。
    """
    segs_gray: List[List[List[float]]] = []
    segs_black: List[List[List[float]]] = []

    def add_seg(dst: List, p0: Tuple[float, float, float], p1: Tuple[float, float, float]) -> None:
        dst.append([list(p0), list(p1)])

    frame_height = 620.0
    frame_half_height = frame_height / 2.0
    frame_face_to_face = 622.0
    frame_radius = frame_face_to_face / 2.0
    n_oct = 8
    oct_vertex_r = frame_radius / math.cos(math.pi / 8.0)
    oct_x_o: List[float] = []
    oct_z_o: List[float] = []
    for i in range(n_oct):
        ang = i * math.pi / 4.0 + math.pi / 8.0
        oct_x_o.append(oct_vertex_r * math.cos(ang))
        oct_z_o.append(oct_vertex_r * math.sin(ang))

    for y_sign in (-1, 1):
        y_o = y_sign * frame_half_height
        for i in range(n_oct):
            j = (i + 1) % n_oct
            xa, ya, za = oct_x_o[i], y_o, oct_z_o[i]
            xb, yb, zb = oct_x_o[j], y_o, oct_z_o[j]
            p0 = tpc_local_to_display(xa, ya, za)
            p1 = tpc_local_to_display(xb, yb, zb)
            add_seg(segs_gray, p0, p1)

    for i in range(n_oct):
        xb, zb = oct_x_o[i], oct_z_o[i]
        p0 = tpc_local_to_display(xb, -frame_half_height, zb)
        p1 = tpc_local_to_display(xb, frame_half_height, zb)
        add_seg(segs_gray, p0, p1)

    target_z_o = -143.0
    target_h = 100.0
    target_half_h = target_h / 2.0
    target_r = 40.0
    n_circle = 64
    n_vert = 8

    for y_sign in (-1, 1):
        y_o = y_sign * target_half_h
        for i in range(n_circle):
            j = (i + 1) % n_circle
            t0 = 2.0 * math.pi * i / n_circle
            t1 = 2.0 * math.pi * j / n_circle
            x0 = target_r * math.cos(t0)
            z0 = target_r * math.sin(t0) + target_z_o
            x1 = target_r * math.cos(t1)
            z1 = target_r * math.sin(t1) + target_z_o
            p0 = tpc_local_to_display(x0, y_o, z0)
            p1 = tpc_local_to_display(x1, y_o, z1)
            add_seg(segs_black, p0, p1)

    for i in range(n_vert):
        th = 2.0 * math.pi * i / n_vert
        x0 = target_r * math.cos(th)
        z0 = target_r * math.sin(th) + target_z_o
        p0 = tpc_local_to_display(x0, -target_half_h, z0)
        p1 = tpc_local_to_display(x0, target_half_h, z0)
        add_seg(segs_black, p0, p1)

    return np.asarray(segs_gray, dtype=float), np.asarray(segs_black, dtype=float)


_TPC_FRAME_SEGMENTS_GRAY, _TPC_FRAME_SEGMENTS_BLACK = _make_tpc_frame_segment_arrays()


def _draw_tpc_frame(ax) -> None:
    """check_tpc_track_3d.C と同様の八角柱＋ターゲット円筒。"""
    if len(_TPC_FRAME_SEGMENTS_GRAY):
        ax.add_collection3d(
            Line3DCollection(
                _TPC_FRAME_SEGMENTS_GRAY,
                colors="0.45",
                linewidths=1.8,
                alpha=0.75,
            )
        )
    if len(_TPC_FRAME_SEGMENTS_BLACK):
        ax.add_collection3d(
            Line3DCollection(
                _TPC_FRAME_SEGMENTS_BLACK,
                colors="k",
                linewidths=1.8,
                alpha=0.9,
            )
        )


def _merge_span_with_frame(
    xmin: float, xmax: float, ymin: float, ymax: float, zmin: float, zmax: float
) -> Tuple[float, float, float, float, float, float]:
    """軸レンジに check_tpc_track_3d.C と同じ表示ボックスを含める。"""
    xmin = min(xmin, DISPLAY_XLIM[0])
    xmax = max(xmax, DISPLAY_XLIM[1])
    ymin = min(ymin, DISPLAY_YLIM[0])
    ymax = max(ymax, DISPLAY_YLIM[1])
    zmin = min(zmin, DISPLAY_ZLIM[0])
    zmax = max(zmax, DISPLAY_ZLIM[1])
    return xmin, xmax, ymin, ymax, zmin, zmax


def helix_theta_draw_range(helix_t_row: List[float]) -> Optional[Tuple[float, float]]:
    """
    ``check_tpc_helix_track_3d.C`` の ``HelixThetaDrawRange`` と同じロジック
    （先頭末尾優先の拡張・パーセンタイル間引き・端マージン）。
    """
    finite: List[float] = []
    gmin = math.inf
    gmax = -math.inf
    t_first: Optional[float] = None
    t_last: Optional[float] = None
    have_first = False
    for w in helix_t_row:
        wf = float(w)
        if not math.isfinite(wf):
            continue
        finite.append(wf)
        if not have_first:
            t_first = wf
            have_first = True
        t_last = wf
        gmin = min(gmin, wf)
        gmax = max(gmax, wf)
    if not have_first or not finite or t_first is None or t_last is None:
        return None
    gspan = gmax - gmin
    if gspan < 1e-12:
        return None
    e_lo = min(t_first, t_last)
    e_hi = max(t_first, t_last)
    espan = e_hi - e_lo
    two_pi = 2.0 * math.pi
    prefer_wide = gspan > two_pi and espan + 1e-9 < gspan * 0.75 and espan > 1e-9
    prefer_out = espan > 1e-9 and espan + 1e-9 < gspan * 0.9 and gspan > espan * 2.5 + 0.3
    if prefer_wide or prefer_out:
        t_lo, t_hi = e_lo, e_hi
    else:
        t_lo, t_hi = gmin, gmax

    if len(finite) >= 5:
        sv = sorted(finite)
        n = len(sv)
        k0 = int(0.06 * (n - 1))
        k1 = int(0.94 * (n - 1))
        p_lo, p_hi = sv[k0], sv[k1]
        pspan = p_hi - p_lo
        cspan = t_hi - t_lo
        if pspan > 1e-9 and cspan > pspan * 1.3 + 0.2:
            t_lo, t_hi = p_lo, p_hi

    mrg = (t_hi - t_lo) * 0.02 + 5e-4
    t_lo -= mrg
    t_hi += mrg
    if t_hi - t_lo < 1e-12:
        return None
    return t_lo, t_hi


def helix_polyline_sample_count(t0: float, t1: float, cap: int = 512) -> int:
    """θ 区間が長いほどポリライン点数を増やす（短い区間で棒状に見えるのを防ぐ）。"""
    span = abs(t1 - t0)
    n = int(48.0 + 64.0 * span / (2.0 * math.pi))
    return max(72, min(cap, n))


def helix_draw_theta_range(ev: EventDataHelix3D, itrack: int) -> Optional[Tuple[float, float]]:
    """
    螺旋ポリライン用の (t_lo, t_hi)。枝 ``helix_theta_min`` / ``helix_theta_max``（または旧名）があれば
    GetMint/GetMaxt に相当する範囲を優先（表示用にわずかにマージン）。
    無い場合は ``helix_t`` から ``helix_theta_draw_range`` で推定。
    """
    if (
        g_has_helix_t_range
        and itrack < len(ev.helix_t_min)
        and itrack < len(ev.helix_t_max)
    ):
        a = float(ev.helix_t_min[itrack])
        b = float(ev.helix_t_max[itrack])
        if math.isfinite(a) and math.isfinite(b) and abs(b - a) > 1e-12:
            lo, hi = min(a, b), max(a, b)
            mrg = (hi - lo) * 0.02 + 5e-4
            return lo - mrg, hi + mrg
    if itrack < len(ev.helix_t):
        return helix_theta_draw_range(ev.helix_t[itrack])
    return None


def _set_path_uproot(path: str) -> bool:
    global g_uproot_file, g_uproot_tree, g_entries_helix, g_backend
    global g_has_calpos, g_has_vertex, g_has_close_dist, g_branch_isbeam, g_has_helix_t_range, g_helix_trange_src

    import uproot

    g_uproot_file = uproot.open(path)
    g_uproot_tree = None
    try:
        g_uproot_tree = g_uproot_file["tpc"]
    except KeyError:
        for kk in g_uproot_file.keys():
            base = str(kk).split("/")[-1].split(";")[0]
            if base == "tpc":
                g_uproot_tree = g_uproot_file[kk]
                break
    if g_uproot_tree is None:
        print("Error: Cannot find tree 'tpc' in file")
        try:
            g_uproot_file.close()
        except Exception:
            pass
        g_uproot_file = None
        g_has_helix_t_range = False
        g_helix_trange_src = 0
        return False
    g_entries_helix = int(getattr(g_uproot_tree, "num_entries", len(g_uproot_tree)))
    g_backend = "uproot"

    t = g_uproot_tree
    g_has_calpos = _uproot_has_branch(t, "calpos_x")
    g_has_vertex = _uproot_has_branch(t, "vtxTpc")
    g_has_close_dist = _uproot_has_branch(t, "closeDistTpc")
    g_branch_isbeam = _uproot_has_branch(t, "isBeam")
    if _uproot_has_branch(t, "helix_theta_min") and _uproot_has_branch(t, "helix_theta_max"):
        g_helix_trange_src = 1
    elif _uproot_has_branch(t, "helix_t_min") and _uproot_has_branch(t, "helix_t_max"):
        g_helix_trange_src = 2
    else:
        g_helix_trange_src = 0
    g_has_helix_t_range = g_helix_trange_src != 0

    required = (
        "ntTpc",
        "nhtrack",
        "helix_cx",
        "helix_cy",
        "helix_z0",
        "helix_r",
        "helix_dz",
        "helix_t",
        "hitpos_x",
        "hitpos_y",
        "hitpos_z",
        "pid",
    )
    missing = [b for b in required if not _uproot_has_branch(t, b)]
    if missing:
        print(f"Error: tree 'tpc' is missing required branch(es): {', '.join(missing)}")
        print("       (This file may not be from DstTPCHelixTracking.)")
        try:
            g_uproot_file.close()
        except Exception:
            pass
        g_uproot_file = None
        g_uproot_tree = None
        g_entries_helix = 0
        g_backend = ""
        g_has_helix_t_range = False
        g_helix_trange_src = 0
        return False
    return True


def _set_path_pyroot(path: str) -> bool:
    global g_file_helix, g_tree_helix, g_entries_helix, g_backend
    global g_has_calpos, g_has_vertex, g_has_close_dist, g_branch_isbeam, g_has_helix_t_range, g_helix_trange_src

    ROOT = _ensure_pyroot()
    g_file_helix = ROOT.TFile.Open(path)
    if not g_file_helix or g_file_helix.IsZombie():
        print(f"Error: Cannot open file {path}")
        return False

    g_tree_helix = g_file_helix.Get("tpc")
    if not g_tree_helix:
        print("Error: Cannot find tree 'tpc' in file")
        return False

    g_entries_helix = int(g_tree_helix.GetEntries())
    g_backend = "pyroot"

    g_has_calpos = bool(g_tree_helix.GetBranch("calpos_x"))
    g_has_vertex = bool(g_tree_helix.GetBranch("vtxTpc"))
    g_has_close_dist = bool(g_tree_helix.GetBranch("closeDistTpc"))
    g_branch_isbeam = bool(g_tree_helix.GetBranch("isBeam"))
    if g_tree_helix.GetBranch("helix_theta_min") and g_tree_helix.GetBranch("helix_theta_max"):
        g_helix_trange_src = 1
    elif g_tree_helix.GetBranch("helix_t_min") and g_tree_helix.GetBranch("helix_t_max"):
        g_helix_trange_src = 2
    else:
        g_helix_trange_src = 0
    g_has_helix_t_range = g_helix_trange_src != 0

    required = (
        "ntTpc",
        "nhtrack",
        "helix_cx",
        "helix_cy",
        "helix_z0",
        "helix_r",
        "helix_dz",
        "helix_t",
        "hitpos_x",
        "hitpos_y",
        "hitpos_z",
        "pid",
    )
    missing = [b for b in required if not g_tree_helix.GetBranch(b)]
    if missing:
        print(f"Error: tree 'tpc' is missing required branch(es): {', '.join(missing)}")
        print("       (This file may not be from DstTPCHelixTracking.)")
        try:
            g_file_helix.Close()
        except Exception:
            pass
        g_file_helix = None
        g_tree_helix = None
        g_entries_helix = 0
        g_backend = ""
        g_has_helix_t_range = False
        g_helix_trange_src = 0
        return False
    return True


def set_path_helix_3d(path: str, read_isbeam: bool = False, backend: str = "auto") -> None:
    global g_read_isbeam, g_branch_isbeam

    _close_helix_sources()

    g_read_isbeam = bool(read_isbeam)
    be = (backend or "auto").strip().lower()
    if be not in ("auto", "uproot", "pyroot"):
        print(f"Error: unknown backend {backend!r} (use auto, uproot, pyroot)")
        return

    ok = False
    if be == "pyroot":
        ok = _set_path_pyroot(path)
    elif be == "uproot":
        try:
            ok = _set_path_uproot(path)
        except ImportError:
            print("Error: uproot/awkward が入っていません: pip install uproot awkward")
            ok = False
    else:
        try:
            ok = _set_path_uproot(path)
        except ImportError:
            print("Info: uproot 未使用 → PyROOT で開きます（セグフォする場合は pip install uproot awkward）")
            ok = _set_path_pyroot(path)

    if not ok:
        g_read_isbeam = False
        return

    if g_read_isbeam and not g_branch_isbeam:
        print("Warning: read_isbeam=True but tree has no 'isBeam' branch; treating as absent.")
        g_read_isbeam = False

    print(f"File opened: {path}  [backend={g_backend}]")
    print(f"Total entries: {g_entries_helix}")
    print(
        f"Optional branches: calpos={g_has_calpos}, vtx={g_has_vertex}, "
        f"closeDist={g_has_close_dist}, helix_theta_range={g_has_helix_t_range}(mode={g_helix_trange_src}), "
        f"isBeam_read={g_read_isbeam and g_branch_isbeam}"
    )


def _load_event_helix_uproot(entry: int) -> Optional[int]:
    global g_event_helix

    t = g_uproot_tree
    base = [
        "run_number",
        "event_number",
        "ntTpc",
        "nhtrack",
        "pid",
        "helix_cx",
        "helix_cy",
        "helix_z0",
        "helix_r",
        "helix_dz",
        "helix_t",
        "hitpos_x",
        "hitpos_y",
        "hitpos_z",
    ]
    extra: List[str] = []
    if g_has_calpos:
        extra += ["calpos_x", "calpos_y", "calpos_z"]
    if g_has_vertex:
        extra += ["vtxTpc", "vtyTpc", "vtzTpc"]
    if g_has_close_dist:
        extra.append("closeDistTpc")
    if g_read_isbeam and g_branch_isbeam:
        extra.append("isBeam")
    if _uproot_has_branch(t, "charge"):
        extra.append("charge")
    if _uproot_has_branch(t, "mom0"):
        extra.append("mom0")
    if g_helix_trange_src == 1:
        extra += ["helix_theta_min", "helix_theta_max"]
    elif g_helix_trange_src == 2:
        extra += ["helix_t_min", "helix_t_max"]

    names = [n for n in base + extra if _uproot_has_branch(t, n)]
    try:
        chunk = t.arrays(names, library="ak", entry_start=entry, entry_stop=entry + 1)
    except (ImportError, ModuleNotFoundError):
        print("Error: uproot でネスト枝を読むには awkward が必要です: pip install awkward")
        return None

    ev = g_event_helix
    ev.runnum = _cell0_scalar_int(chunk, "run_number") if "run_number" in names else 0
    ev.evnum = _cell0_scalar_int(chunk, "event_number") if "event_number" in names else 0
    ev.ntTpc = _cell0_scalar_int(chunk, "ntTpc")
    nt = ev.ntTpc

    ev.nhtrack = _cell0_vec_int(chunk, "nhtrack")
    ev.pid = _cell0_vec_int(chunk, "pid")
    ev.is_beam = _cell0_vec_int(chunk, "isBeam") if (g_read_isbeam and g_branch_isbeam and "isBeam" in names) else []
    ev.charge = _cell0_vec_int(chunk, "charge") if "charge" in names else []
    ev.mom0 = _cell0_vec_float(chunk, "mom0") if "mom0" in names else []

    ev.helix_cx = _cell0_vec_float(chunk, "helix_cx")
    ev.helix_cy = _cell0_vec_float(chunk, "helix_cy")
    ev.helix_z0 = _cell0_vec_float(chunk, "helix_z0")
    ev.helix_r = _cell0_vec_float(chunk, "helix_r")
    ev.helix_dz = _cell0_vec_float(chunk, "helix_dz")
    if g_helix_trange_src == 1 and "helix_theta_min" in names and "helix_theta_max" in names:
        ev.helix_t_min = _cell0_vec_float(chunk, "helix_theta_min")
        ev.helix_t_max = _cell0_vec_float(chunk, "helix_theta_max")
    elif g_helix_trange_src == 2 and "helix_t_min" in names and "helix_t_max" in names:
        ev.helix_t_min = _cell0_vec_float(chunk, "helix_t_min")
        ev.helix_t_max = _cell0_vec_float(chunk, "helix_t_max")
    else:
        ev.helix_t_min = []
        ev.helix_t_max = []

    ev.helix_t = _cell0_nested_float(chunk, "helix_t", nt)
    ev.hitpos_x = _cell0_nested_float(chunk, "hitpos_x", nt)
    ev.hitpos_y = _cell0_nested_float(chunk, "hitpos_y", nt)
    ev.hitpos_z = _cell0_nested_float(chunk, "hitpos_z", nt)

    ev.calpos_x = []
    ev.calpos_y = []
    ev.calpos_z = []
    if g_has_calpos and "calpos_x" in names:
        ev.calpos_x = _cell0_nested_float(chunk, "calpos_x", nt)
        ev.calpos_y = _cell0_nested_float(chunk, "calpos_y", nt)
        ev.calpos_z = _cell0_nested_float(chunk, "calpos_z", nt)

    ev.vtxTpc = []
    ev.vtyTpc = []
    ev.vtzTpc = []
    ev.closeDistTpc = []
    if g_has_vertex and "vtxTpc" in names:
        ev.vtxTpc = _cell0_nested_float(chunk, "vtxTpc", nt)
        ev.vtyTpc = _cell0_nested_float(chunk, "vtyTpc", nt)
        ev.vtzTpc = _cell0_nested_float(chunk, "vtzTpc", nt)
    if g_has_close_dist and "closeDistTpc" in names:
        ev.closeDistTpc = _cell0_nested_float(chunk, "closeDistTpc", nt)

    print(f"Event loaded: Run={ev.runnum}, Event={ev.evnum}, Entry={entry}")
    print(f"  Tracks: {ev.ntTpc}")
    return entry


def _load_event_helix_pyroot(entry: int) -> Optional[int]:
    global g_tree_helix, g_event_helix

    g_tree_helix.GetEntry(entry)

    ev = g_event_helix
    ev.runnum = int(getattr(g_tree_helix, "run_number", 0))
    ev.evnum = int(getattr(g_tree_helix, "event_number", 0))
    ev.ntTpc = int(g_tree_helix.ntTpc)
    nt = ev.ntTpc

    ev.nhtrack = _vec_int_to_list(g_tree_helix.nhtrack)
    if g_read_isbeam and g_branch_isbeam:
        ev.is_beam = _vec_int_to_list(g_tree_helix.isBeam)
    else:
        ev.is_beam = []
    ev.charge = _vec_int_to_list(g_tree_helix.charge) if g_tree_helix.GetBranch("charge") else []
    ev.mom0 = _vec_double_to_list(g_tree_helix.mom0) if g_tree_helix.GetBranch("mom0") else []
    ev.pid = _vec_int_to_list(g_tree_helix.pid)

    ev.helix_cx = _vec_double_to_list(g_tree_helix.helix_cx)
    ev.helix_cy = _vec_double_to_list(g_tree_helix.helix_cy)
    ev.helix_z0 = _vec_double_to_list(g_tree_helix.helix_z0)
    ev.helix_r = _vec_double_to_list(g_tree_helix.helix_r)
    ev.helix_dz = _vec_double_to_list(g_tree_helix.helix_dz)
    if g_helix_trange_src == 1:
        ev.helix_t_min = _vec_double_to_list(g_tree_helix.helix_theta_min)
        ev.helix_t_max = _vec_double_to_list(g_tree_helix.helix_theta_max)
    elif g_helix_trange_src == 2:
        ev.helix_t_min = _vec_double_to_list(g_tree_helix.helix_t_min)
        ev.helix_t_max = _vec_double_to_list(g_tree_helix.helix_t_max)
    else:
        ev.helix_t_min = []
        ev.helix_t_max = []

    ev.helix_t = _nested_vec_double_from_tree(g_tree_helix, "helix_t", nt)
    ev.hitpos_x = _nested_vec_double_from_tree(g_tree_helix, "hitpos_x", nt)
    ev.hitpos_y = _nested_vec_double_from_tree(g_tree_helix, "hitpos_y", nt)
    ev.hitpos_z = _nested_vec_double_from_tree(g_tree_helix, "hitpos_z", nt)

    ev.calpos_x = []
    ev.calpos_y = []
    ev.calpos_z = []
    if g_has_calpos:
        ev.calpos_x = _nested_vec_double_from_tree(g_tree_helix, "calpos_x", nt)
        ev.calpos_y = _nested_vec_double_from_tree(g_tree_helix, "calpos_y", nt)
        ev.calpos_z = _nested_vec_double_from_tree(g_tree_helix, "calpos_z", nt)

    ev.vtxTpc = []
    ev.vtyTpc = []
    ev.vtzTpc = []
    ev.closeDistTpc = []
    if g_has_vertex:
        ev.vtxTpc = _nested_vec_double_from_tree(g_tree_helix, "vtxTpc", nt)
        ev.vtyTpc = _nested_vec_double_from_tree(g_tree_helix, "vtyTpc", nt)
        ev.vtzTpc = _nested_vec_double_from_tree(g_tree_helix, "vtzTpc", nt)
    if g_has_close_dist:
        ev.closeDistTpc = _nested_vec_double_from_tree(g_tree_helix, "closeDistTpc", nt)

    print(f"Event loaded: Run={ev.runnum}, Event={ev.evnum}, Entry={entry}")
    print(f"  Tracks: {ev.ntTpc}")
    return entry


def load_event_helix_3d(entry: int = -1) -> Optional[int]:
    """
    イベントを読み込む。成功時は実際に読んだツリーエントリ番号を返す。失敗時は None。
    entry < 0 のときはランダムに選ぶ。
    """
    global g_entries_helix, g_event_helix, g_current_helix_entry

    if not _io_ready():
        print("Error: No file opened. Use set_path_helix_3d() first.")
        return None
    if g_entries_helix == 0:
        print("Error: No entries in tree")
        return None

    if entry < 0:
        entry = random.randint(0, g_entries_helix - 1)
    if entry >= g_entries_helix:
        print(f"Error: Entry {entry} is out of range [0, {g_entries_helix})")
        return None

    if g_backend == "uproot":
        res = _load_event_helix_uproot(entry)
    else:
        res = _load_event_helix_pyroot(entry)
    if res is not None:
        g_current_helix_entry = res
    return res


def find_entry_helix_3d(
    run: int,
    event: Optional[int] = None,
    *,
    start_entry: int = 0,
    max_scan: Optional[int] = None,
) -> Optional[int]:
    """
    ``run_number`` が ``run`` に一致する最初のツリーエントリを返す。
    ``event`` を与えたときは ``event_number`` も一致が必要。

    ``max_scan`` を正の整数にすると ``start_entry`` からその件数だけを走査する
    （見つからなければ None）。None または 0 以下はファイル末尾まで走査。
    """
    if not _io_ready():
        print("Error: No file opened. Use set_path_helix_3d() first.")
        return None
    if start_entry < 0 or start_entry >= g_entries_helix:
        print(f"Error: start_entry {start_entry} out of range [0, {g_entries_helix})")
        return None

    scan_end = g_entries_helix
    if max_scan is not None and max_scan > 0:
        scan_end = min(start_entry + max_scan, g_entries_helix)

    if g_backend == "uproot":
        t = g_uproot_tree
        if not (_uproot_has_branch(t, "run_number") and _uproot_has_branch(t, "event_number")):
            print("Error: tree has no run_number / event_number for scan.")
            return None
        batch = 65536
        es = start_entry
        while es < scan_end:
            ee = min(es + batch, scan_end)
            try:
                chunk = t.arrays(
                    ["run_number", "event_number"],
                    entry_start=es,
                    entry_stop=ee,
                    library="np",
                )
            except Exception as ex:
                print(f"Error reading run/event columns: {ex}")
                return None
            runs = np.asarray(chunk["run_number"], dtype=np.int64)
            evs = np.asarray(chunk["event_number"], dtype=np.int64)
            mask = runs == int(run)
            if event is not None:
                mask &= evs == int(event)
            hit = np.flatnonzero(mask)
            if hit.size > 0:
                found = es + int(hit[0])
                print(f"find_entry_helix_3d: run={run} event={event} -> entry={found}")
                return found
            es = ee
    else:
        tr = g_tree_helix
        if tr is None:
            return None
        if not tr.GetBranch("run_number") or not tr.GetBranch("event_number"):
            print("Error: tree has no run_number / event_number for scan.")
            return None
        for e in range(start_entry, scan_end):
            tr.GetEntry(e)
            rn = int(getattr(tr, "run_number", -1))
            if rn != int(run):
                continue
            en = int(getattr(tr, "event_number", -1))
            if event is not None and en != int(event):
                continue
            print(f"find_entry_helix_3d: run={run} event={event} -> entry={e}")
            return e

    print(
        f"find_entry_helix_3d: no match for run={run} event={event} "
        f"in entries [{start_entry}, {scan_end})"
    )
    return None


def find_entry_by_event_number(
    event_num: int,
    *,
    start_entry: int = 0,
    max_scan: Optional[int] = None,
) -> Optional[int]:
    """
    ``event_number`` が ``event_num`` に一致する最初のツリーエントリを返す（run は見ない）。

    複数ランで同じ event_number が繰り返される場合は **先頭に近い方** が選ばれる。
    ``max_scan`` で走査長を制限可能（``find_entry_helix_3d`` と同様）。
    """
    if not _io_ready():
        print("Error: No file opened. Use set_path_helix_3d() first.")
        return None
    if start_entry < 0 or start_entry >= g_entries_helix:
        print(f"Error: start_entry {start_entry} out of range [0, {g_entries_helix})")
        return None

    scan_end = g_entries_helix
    if max_scan is not None and max_scan > 0:
        scan_end = min(start_entry + max_scan, g_entries_helix)

    if g_backend == "uproot":
        t = g_uproot_tree
        if not _uproot_has_branch(t, "event_number"):
            print("Error: tree has no event_number for scan.")
            return None
        batch = 65536
        es = start_entry
        want = int(event_num)
        while es < scan_end:
            ee = min(es + batch, scan_end)
            try:
                chunk = t.arrays(
                    ["event_number"],
                    entry_start=es,
                    entry_stop=ee,
                    library="np",
                )
            except Exception as ex:
                print(f"Error reading event_number column: {ex}")
                return None
            evs = np.asarray(chunk["event_number"], dtype=np.int64)
            hit = np.flatnonzero(evs == want)
            if hit.size > 0:
                found = es + int(hit[0])
                print(f"find_entry_by_event_number: event_number={event_num} -> entry={found}")
                return found
            es = ee
    else:
        tr = g_tree_helix
        if tr is None:
            return None
        if not tr.GetBranch("event_number"):
            print("Error: tree has no event_number for scan.")
            return None
        want = int(event_num)
        for e in range(start_entry, scan_end):
            tr.GetEntry(e)
            if int(getattr(tr, "event_number", -1)) == want:
                print(f"find_entry_by_event_number: event_number={event_num} -> entry={e}")
                return e

    print(
        f"find_entry_by_event_number: no match for event_number={event_num} "
        f"in entries [{start_entry}, {scan_end})"
    )
    return None


def _load_scan_only_uproot(entry: int) -> bool:
    global g_event_helix, g_uproot_tree

    t = g_uproot_tree
    cols = ["ntTpc"]
    if g_has_close_dist:
        cols.append("closeDistTpc")
    if g_read_isbeam and g_branch_isbeam:
        cols.append("isBeam")
    cols = [c for c in cols if _uproot_has_branch(t, c)]
    try:
        chunk = t.arrays(cols, library="ak", entry_start=entry, entry_stop=entry + 1)
    except (ImportError, ModuleNotFoundError):
        print("Error: awkward が必要です: pip install awkward")
        return False
    ev = g_event_helix
    ev.ntTpc = _cell0_scalar_int(chunk, "ntTpc")
    ev.is_beam = _cell0_vec_int(chunk, "isBeam") if "isBeam" in cols else []
    ev.closeDistTpc = _cell0_nested_float(chunk, "closeDistTpc", ev.ntTpc) if "closeDistTpc" in cols else []
    return True


def _load_scan_only_pyroot(entry: int) -> bool:
    global g_tree_helix, g_event_helix

    if not g_tree_helix:
        return False
    if entry < 0 or entry >= g_entries_helix:
        return False

    g_tree_helix.GetEntry(entry)
    ev = g_event_helix
    ev.ntTpc = int(g_tree_helix.ntTpc)
    ev.is_beam = []
    if g_read_isbeam and g_branch_isbeam:
        ev.is_beam = _vec_int_to_list(g_tree_helix.isBeam)
    ev.closeDistTpc = []
    if g_has_close_dist:
        ev.closeDistTpc = _nested_vec_double_from_tree(g_tree_helix, "closeDistTpc", ev.ntTpc)
    return True


def _load_scan_only(entry: int) -> bool:
    """scan_kp_candidate 専用: ntTpc / isBeam(任意) / closeDistTpc だけを読む。"""
    if not _io_ready() or entry < 0 or entry >= g_entries_helix:
        return False
    if g_backend == "uproot":
        return _load_scan_only_uproot(entry)
    return _load_scan_only_pyroot(entry)


def _kp_topology_ok(ev: EventDataHelix3D) -> bool:
    if ev.ntTpc != 2:
        return False
    if g_read_isbeam and g_branch_isbeam:
        if len(ev.is_beam) < 2:
            return False
        beams = list(ev.is_beam[:2])
        return beams.count(1) == 1 and beams.count(0) == 1
    # isBeam を読まないときはここでは判定しない（scan 側で closeDist 必須にする）
    return True


def _close_dist_ok(ev: EventDataHelix3D, close_dist_max: Optional[float]) -> bool:
    if close_dist_max is None:
        return True
    if not g_has_close_dist or len(ev.closeDistTpc) < 2:
        return True
    d01 = ev.closeDistTpc[0][1] if len(ev.closeDistTpc[0]) > 1 else float("nan")
    d10 = ev.closeDistTpc[1][0] if len(ev.closeDistTpc[1]) > 0 else float("nan")
    cand = []
    if math.isfinite(d01):
        cand.append(abs(d01))
    if math.isfinite(d10):
        cand.append(abs(d10))
    if not cand:
        return False
    return min(cand) <= close_dist_max


def scan_kp_candidate(
    start_entry: int = 0,
    max_scan: int = 10000,
    close_dist_max: Optional[float] = None,
) -> Optional[int]:
    """
    Kp 弾性散乱の目視用のエントリ走査。

    - ``read_isbeam=True`` で開いていて ``isBeam`` 枝がある場合: 2 トラックかつ
      isBeam が 1 と 0 の組み合わせを要求する。
    - 既定（isBeam を読まない）の場合: 2 トラックに加え、**close_dist_max を必須**
      （最接近距離で散乱らしいイベントに絞る）。isBeam 無しで全 2 トラックを
      走査すると件数が多すぎるため。
    """
    if not _io_ready():
        print("Error: No tree loaded.")
        return None

    use_isbeam_filter = g_read_isbeam and g_branch_isbeam
    if not use_isbeam_filter and close_dist_max is None:
        print(
            "scan_kp_candidate: isBeam is not read (default, avoids segfault on broken branch). "
            "Either open with read_isbeam=True, or pass close_dist_max=... to scan 2-track+vertex distance."
        )
        return None

    end = min(start_entry + max_scan, g_entries_helix)
    for e in range(start_entry, end):
        if not _load_scan_only(e):
            continue
        ev = g_event_helix
        if not _kp_topology_ok(ev):
            continue
        if not _close_dist_ok(ev, close_dist_max):
            continue
        print(f"scan_kp_candidate: found entry {e} (scanned {e - start_entry + 1} events)")
        return e
    print(f"scan_kp_candidate: no candidate in [{start_entry}, {end})")
    return None


def _collect_bounds(ev: EventDataHelix3D) -> Tuple[float, float, float, float, float, float]:
    """表示座標 (z,x,y) でヒット＋螺旋のバウンディング。"""
    xs, ys, zs = [], [], []
    for i in range(ev.ntTpc):
        pid_i = ev.pid[i] if i < len(ev.pid) else 0
        if not track_visible_kp_filter(pid_i):
            continue
        for j in range(len(ev.hitpos_x[i])):
            mx, my, mz = tpc_local_to_display(ev.hitpos_x[i][j], ev.hitpos_y[i][j], ev.hitpos_z[i][j])
            xs.append(mx)
            ys.append(my)
            zs.append(mz)
        tr = helix_draw_theta_range(ev, i)
        rr_b = ev.helix_r[i] if i < len(ev.helix_r) else 0.0
        if tr and i < len(ev.helix_cx) and math.isfinite(float(rr_b)):
            t0, t1 = tr
            nt_b = min(200, max(48, helix_polyline_sample_count(t0, t1, cap=200)))
            th = np.linspace(t0, t1, nt_b)
            x, y, z = helix_xyz(ev.helix_cx[i], ev.helix_cy[i], ev.helix_z0[i], ev.helix_r[i], ev.helix_dz[i], th)
            mx, my, mz = tpc_local_to_display_vec(x, y, z)
            xs.extend(mx.tolist())
            ys.extend(my.tolist())
            zs.extend(mz.tolist())
    if not xs:
        return (
            DISPLAY_XLIM[0],
            DISPLAY_XLIM[1],
            DISPLAY_YLIM[0],
            DISPLAY_YLIM[1],
            DISPLAY_ZLIM[0],
            DISPLAY_ZLIM[1],
        )
    margin = 0.1

    def span(a: List[float]) -> Tuple[float, float]:
        lo, hi = min(a), max(a)
        m = (hi - lo) * margin if hi > lo else 10.0
        return lo - m, hi + m

    xmin, xmax = span(xs)
    ymin, ymax = span(ys)
    zmin, zmax = span(zs)
    return xmin, xmax, ymin, ymax, zmin, zmax


def _style_ax_without_builtin_axes(ax) -> None:
    """既定の 3D 軸パネル・目盛・軸ラインを非表示（座標ガイドは手描き）。"""
    ax.set_axis_off()
    try:
        ax.grid(False)
    except Exception:
        pass


def _draw_manual_xyz_axes(
    ax,
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    zmin: float,
    zmax: float,
) -> None:
    """
    check_tpc_track_3d.C と同様の角ガイド。表示座標 (Xv,Yv,Zv)=(局所z,局所x,局所y) のため、
    ラベルは局所軸名 z, x, y（色は赤・緑・青で各矢印に対応）。
    """
    sx = xmax - xmin
    sy = ymax - ymin
    sz = zmax - zmin
    ox = xmin + 0.05 * sx
    oy = ymin + 0.05 * sy
    oz = zmin + 0.05 * sz
    L = min(sx, sy, sz) * 0.1
    if L <= 0:
        L = 20.0
    arr = L * 0.15
    lbl_off = L * 0.22
    fs = 9

    ax.plot([ox, ox + L], [oy, oy], [oz, oz], color="tab:red", linewidth=2.0, zorder=8)
    ax.plot(
        [ox + L, ox + L - arr],
        [oy, oy - arr * 0.3],
        [oz, oz],
        color="tab:red",
        linewidth=2.0,
        zorder=8,
    )
    ax.plot(
        [ox + L, ox + L - arr],
        [oy, oy + arr * 0.3],
        [oz, oz],
        color="tab:red",
        linewidth=2.0,
        zorder=8,
    )
    ax.text(ox + L + lbl_off, oy, oz, "z", color="tab:red", fontsize=fs, va="center", ha="left")

    ax.plot([ox, ox], [oy, oy + L], [oz, oz], color="tab:green", linewidth=2.0, zorder=8)
    ax.plot(
        [ox, ox - arr * 0.3],
        [oy + L, oy + L - arr],
        [oz, oz],
        color="tab:green",
        linewidth=2.0,
        zorder=8,
    )
    ax.plot(
        [ox, ox + arr * 0.3],
        [oy + L, oy + L - arr],
        [oz, oz],
        color="tab:green",
        linewidth=2.0,
        zorder=8,
    )
    ax.text(ox, oy + L + lbl_off, oz, "x", color="tab:green", fontsize=fs, va="bottom", ha="center")

    ax.plot([ox, ox], [oy, oy], [oz, oz + L], color="tab:blue", linewidth=2.0, zorder=8)
    ax.plot(
        [ox, ox - arr * 0.3],
        [oy, oy],
        [oz + L, oz + L - arr],
        color="tab:blue",
        linewidth=2.0,
        zorder=8,
    )
    ax.plot(
        [ox, ox + arr * 0.3],
        [oy, oy],
        [oz + L, oz + L - arr],
        color="tab:blue",
        linewidth=2.0,
        zorder=8,
    )
    ax.text(ox, oy, oz + L + lbl_off, "y", color="tab:blue", fontsize=fs, va="bottom", ha="center")


def _draw_vertex(ax, ev: EventDataHelix3D) -> None:
    if ev.ntTpc != 2 or not g_has_vertex:
        return
    if g_helix_kp_only:
        p0 = ev.pid[0] if len(ev.pid) > 0 else 0
        p1 = ev.pid[1] if len(ev.pid) > 1 else 0
        if not (track_visible_kp_filter(p0) and track_visible_kp_filter(p1)):
            return
    if len(ev.vtxTpc) < 2:
        return
    vx = ev.vtxTpc[0][1] if len(ev.vtxTpc[0]) > 1 else float("nan")
    vy = ev.vtyTpc[0][1] if len(ev.vtyTpc[0]) > 1 else float("nan")
    vz = ev.vtzTpc[0][1] if len(ev.vtzTpc[0]) > 1 else float("nan")
    if not (math.isfinite(vx) and math.isfinite(vy) and math.isfinite(vz)):
        return
    mx, my, mz = tpc_local_to_display(vx, vy, vz)
    vlabel = None if g_helix_png_export else "TPC vertex (0-1)"
    ax.scatter([mx], [my], [mz], c="magenta", marker="X", s=180, alpha=0.95, label=vlabel)


def helix_default_png_path(ev: EventDataHelix3D) -> str:
    """現在のイベント・エントリから PNG ファイル名を生成する。"""
    ent = g_current_helix_entry if g_current_helix_entry is not None else 0
    return f"helix_run{ev.runnum}_ev{ev.evnum}_entry{ent}.png"


def _helix_apply_transparent_bg(fig, ax) -> None:
    """Figure / 3D Axes の背景・パネルを透過にする（savefig transparent 用）。"""
    fig.patch.set_alpha(0.0)
    fig.patch.set_facecolor("none")
    if ax is None:
        return
    ax.set_facecolor("none")
    try:
        ax.xaxis.pane.set_alpha(0.0)
        ax.yaxis.pane.set_alpha(0.0)
        ax.zaxis.pane.set_alpha(0.0)
        ax.xaxis.pane.set_edgecolor("none")
        ax.yaxis.pane.set_edgecolor("none")
        ax.zaxis.pane.set_edgecolor("none")
    except Exception:
        pass


def _helix_figure_valid() -> bool:
    if g_fig_helix is None:
        return False
    try:
        return bool(plt.fignum_exists(g_fig_helix.number))
    except Exception:
        return False


def _render_helix_3d_plot(ev: EventDataHelix3D, n_helix_samples: int) -> None:
    """現在の g_ax_helix に 1 イベント分を描き直す（plt.show は呼ばない）。"""
    ax = g_ax_helix
    if ax is None:
        return
    ax.clear()

    xmin, xmax, ymin, ymax, zmin, zmax = _collect_bounds(ev)
    xmin, xmax, ymin, ymax, zmin, zmax = _merge_span_with_frame(xmin, xmax, ymin, ymax, zmin, zmax)

    _draw_tpc_frame(ax)

    for itrack in range(ev.ntTpc):
        pid_code = ev.pid[itrack] if itrack < len(ev.pid) else 0
        if not track_visible_kp_filter(pid_code):
            continue
        color, pk = pid_display_color(pid_code)
        lbl = f"helix tr{itrack} {pk} pid=0x{pid_code:x} ({decode_pid_candidates(pid_code)})"
        beam_tag = ""
        if itrack < len(ev.is_beam):
            beam_tag = " beam" if ev.is_beam[itrack] else " scat"
            lbl += beam_tag

        cx = ev.helix_cx[itrack] if itrack < len(ev.helix_cx) else 0.0
        cy = ev.helix_cy[itrack] if itrack < len(ev.helix_cy) else 0.0
        z0 = ev.helix_z0[itrack] if itrack < len(ev.helix_z0) else 0.0
        rr = ev.helix_r[itrack] if itrack < len(ev.helix_r) else 0.0
        dz = ev.helix_dz[itrack] if itrack < len(ev.helix_dz) else 0.0

        tr = helix_draw_theta_range(ev, itrack)
        # 半径が小さくても θ 区間が取れていれば描画（旧 0.05 mm 閾値は線が消える原因のため廃止）
        if tr is not None and math.isfinite(float(rr)):
            t0, t1 = tr
            cap = max(512, n_helix_samples)
            nt = max(8, n_helix_samples, helix_polyline_sample_count(t0, t1, cap=cap))
            if abs(t1 - t0) < 1e-9:
                thetas = np.linspace(t0 - 1e-6, t1 + 1e-6, nt)
            else:
                thetas = np.linspace(t0, t1, nt)
            xh, yh, zh = helix_xyz(cx, cy, z0, rr, dz, thetas)
            mx, my, mz = tpc_local_to_display_vec(xh, yh, zh)
            fin = np.isfinite(mx) & np.isfinite(my) & np.isfinite(mz)
            if np.any(fin):
                ax.plot(
                    mx[fin],
                    my[fin],
                    mz[fin],
                    color=color,
                    linewidth=2.2,
                    linestyle=":",
                    label=None if g_helix_png_export else lbl,
                    zorder=6,
                )

        nh = ev.nhtrack[itrack] if itrack < len(ev.nhtrack) else 0
        if nh > 0 and itrack < len(ev.hitpos_x) and len(ev.hitpos_x[itrack]) > 0:
            nplot = min(nh, len(ev.hitpos_x[itrack]))
            hx = np.asarray(ev.hitpos_x[itrack][:nplot], dtype=float)
            hy = np.asarray(ev.hitpos_y[itrack][:nplot], dtype=float)
            hz = np.asarray(ev.hitpos_z[itrack][:nplot], dtype=float)
            mx, my, mz = tpc_local_to_display_vec(hx, hy, hz)
            ax.scatter(mx, my, mz, color=color, marker="o", s=36, alpha=0.85)

        if g_has_calpos and itrack < len(ev.calpos_x) and len(ev.calpos_x[itrack]) > 0:
            nplot = min(nh, len(ev.calpos_x[itrack])) if nh else len(ev.calpos_x[itrack])
            cxp = np.asarray(ev.calpos_x[itrack][:nplot], dtype=float)
            cyp = np.asarray(ev.calpos_y[itrack][:nplot], dtype=float)
            czp = np.asarray(ev.calpos_z[itrack][:nplot], dtype=float)
            mx, my, mz = tpc_local_to_display_vec(cxp, cyp, czp)
            ax.scatter(mx, my, mz, color=color, marker=".", s=8, alpha=0.35)

    _draw_vertex(ax, ev)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)
    if not g_helix_png_export:
        ent_lbl = g_current_helix_entry if g_current_helix_entry is not None else "?"
        kp_note = "  [--kp: K|p tracks only]" if g_helix_kp_only else ""
        ax.set_title(
            f"TPC Helix 3D (Run {ev.runnum}, Event {ev.evnum}, entry {ent_lbl})  "
            f"[color: K=blue, p=red, pi=green; multi-bit: K>p>pi]  [key n: next entry]{kp_note}"
        )
        ax.legend(loc="upper left", fontsize=8)
    _style_ax_without_builtin_axes(ax)
    _draw_manual_xyz_axes(ax, xmin, xmax, ymin, ymax, zmin, zmax)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    ax.view_init(elev=HELIX_VIEW_ELEV, azim=HELIX_VIEW_AZIM)
    plt.tight_layout()


def _on_helix_key_press(event) -> None:
    global g_current_helix_entry, g_helix_n_samples, g_fig_helix, g_event_helix
    if getattr(event, "key", "") not in ("n", "N"):
        return
    if not _io_ready() or g_current_helix_entry is None:
        return
    if g_fig_helix is None or getattr(event, "canvas", None) is None:
        return
    if event.canvas.figure is not g_fig_helix:
        return
    nxt = g_current_helix_entry + 1
    if nxt >= g_entries_helix:
        nxt = 0
    if load_event_helix_3d(nxt) is None:
        return
    _render_helix_3d_plot(g_event_helix, g_helix_n_samples)
    event.canvas.draw_idle()


def _ensure_helix_figure() -> None:
    global g_fig_helix, g_ax_helix
    if _helix_figure_valid():
        return
    fs = HELIX_SAVE_FIGSIZE_IN if g_helix_png_export else (10.0, 10.0)
    g_fig_helix = plt.figure(figsize=fs)
    g_ax_helix = g_fig_helix.add_subplot(111, projection="3d")
    g_fig_helix.canvas.mpl_connect("key_press_event", _on_helix_key_press)


def event_helix_3d(
    evnum: int = -1,
    n_helix_samples: int = 160,
    *,
    kp_only: bool = False,
    save_png: Optional[str] = None,
    save_dpi: Optional[int] = None,
) -> None:
    global g_helix_n_samples, g_event_helix, g_helix_kp_only, g_helix_png_export, g_fig_helix, g_ax_helix

    g_helix_n_samples = n_helix_samples
    g_helix_kp_only = bool(kp_only)

    if load_event_helix_3d(evnum) is None:
        return

    ev = g_event_helix

    do_save = save_png is not None
    dpi = int(save_dpi) if save_dpi is not None else HELIX_SAVE_DPI
    if dpi < 72:
        dpi = 72

    g_helix_png_export = bool(do_save)
    try:
        _ensure_helix_figure()
        _render_helix_3d_plot(ev, n_helix_samples)

        if do_save:
            sp = (save_png or "").strip()
            out_path = sp if sp else helix_default_png_path(ev)
            out_dir = os.path.dirname(os.path.abspath(out_path))
            if out_dir and not os.path.isdir(out_dir):
                os.makedirs(out_dir, exist_ok=True)
            _helix_apply_transparent_bg(g_fig_helix, g_ax_helix)
            g_fig_helix.savefig(
                out_path,
                dpi=dpi,
                transparent=True,
                facecolor="none",
                edgecolor="none",
                bbox_inches="tight",
                pad_inches=0.06,
                format="png",
            )
            print(f"Saved PNG (transparent, no title/legend): {os.path.abspath(out_path)}  dpi={dpi}", flush=True)
            plt.close(g_fig_helix)
            g_fig_helix = None
            g_ax_helix = None
        else:
            print(
                "  [ヒント] 図ウィンドウをクリックしてフォーカスを当て、キー「n」で次のツリーエントリ",
                flush=True,
            )
            plt.show()
    finally:
        g_helix_png_export = False

    print("\n=== Event Summary (Helix 3D) ===")
    print(f"Tracks: {ev.ntTpc}")
    for itrack in range(ev.ntTpc):
        pid_code = ev.pid[itrack] if itrack < len(ev.pid) else 0
        hid = "" if track_visible_kp_filter(pid_code) else "  [hidden by --kp]"
        ib_str = str(ev.is_beam[itrack]) if itrack < len(ev.is_beam) else "n/a"
        ch = ev.charge[itrack] if itrack < len(ev.charge) else 0
        mom = ev.mom0[itrack] if itrack < len(ev.mom0) else float("nan")
        nh = ev.nhtrack[itrack] if itrack < len(ev.nhtrack) else 0
        print(
            f"  tr{itrack}: isBeam={ib_str} charge={ch} mom0={mom:.4f} nhits={nh} "
            f"pid=0x{pid_code:x} ({decode_pid_candidates(pid_code)}){hid}"
        )
    if ev.ntTpc == 2 and g_has_close_dist and len(ev.closeDistTpc) >= 2:
        d01 = ev.closeDistTpc[0][1] if len(ev.closeDistTpc[0]) > 1 else float("nan")
        print(f"  closeDistTpc[0][1] = {d01}")
    print("================================\n")


def event_helix_3d_browse(first_entry: int = -1, n_helix_samples: int = 160, *, kp_only: bool = False) -> None:
    """
    1 つの図ウィンドウで表示し、キー「n」で次のツリーエントリへ進む。

    ウィンドウを閉じると ``event_helix_3d_browse`` は終了する（従来の「閉じるたび
    に別ウィンドウ」ループは廃止。連続表示は図上の ``n`` を使う）。
    """
    print("Browse: 図をフォーカスしてキー「n」で次エントリ。ウィンドウを閉じて終了。")
    event_helix_3d(first_entry, n_helix_samples, kp_only=kp_only)


def _parse_cli(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="3D event display for DstTPCHelixTracking ROOT output.")
    p.add_argument("rootfile", help="Path to .root file")
    p.add_argument("entry", nargs="?", type=int, default=-1, help="Tree entry (omit or -1 for random)")
    p.add_argument(
        "--scan-kp",
        action="store_true",
        help="Scan for Kp-like events (needs --close-dist-max unless --read-isbeam)",
    )
    p.add_argument("--start", type=int, default=0, help="Start entry for --scan-kp")
    p.add_argument("--max-scan", type=int, default=10000, help="Max entries to scan")
    p.add_argument("--close-dist-max", type=float, default=None, help="Require closeDistTpc <= this (mm)")
    p.add_argument(
        "--read-isbeam",
        action="store_true",
        help="Read isBeam branch (off by default; PyROOT+cppyy では危険なことがある)",
    )
    p.add_argument(
        "--backend",
        choices=("auto", "uproot", "pyroot"),
        default="auto",
        help="auto: uproot+awkward を優先（推奨）。pyroot は --pyroot と同じ。",
    )
    p.add_argument(
        "--pyroot",
        action="store_true",
        help="PyROOT で読む（cppyy セグフォのリスクあり。既定は uproot）",
    )
    p.add_argument(
        "--browse",
        action="store_true",
        help="図ウィンドウでキー「n」で次エントリ。ウィンドウを閉じて終了",
    )
    p.add_argument(
        "--run",
        type=int,
        default=None,
        help="run_number が一致する最初のツリーエントリを開く（任意で --event と併用）",
    )
    p.add_argument(
        "--event",
        type=int,
        default=None,
        help="(1) --run あり: event_number も一致 (2) --run なし: event_number のみで先頭から検索",
    )
    p.add_argument(
        "--max-run-scan",
        type=int,
        default=None,
        help="--run または --event のみ検索で、先頭から走査する最大エントリ数（省略時は末尾まで）",
    )
    p.add_argument(
        "--kp",
        action="store_true",
        help="dE/dx PID で Kaon(0x2) または Proton(0x4) ビットを含むトラックだけ描画",
    )
    p.add_argument(
        "--save",
        nargs="?",
        const="",
        default=None,
        metavar="PATH",
        help="高解像度・透過背景の PNG に保存（タイトル・凡例なし）。PATH 省略時は helix_run*_ev*_entry*.png",
    )
    p.add_argument(
        "--save-dpi",
        type=int,
        default=None,
        metavar="DPI",
        help=f"--save 時の DPI（既定 {HELIX_SAVE_DPI}）",
    )
    return p.parse_args(argv)


def main(argv: List[str]) -> None:
    args = _parse_cli(argv)
    if args.browse and args.save is not None:
        print("注意: --save 指定時は最初の1イベントのみ PNG 保存し、--browse は無視します。", flush=True)
    if args.scan_kp and args.run is not None:
        print("エラー: --scan-kp と --run は同時に使えません。", file=sys.stderr)
        sys.exit(3)
    if args.scan_kp and not args.read_isbeam and args.close_dist_max is None:
        print(
            "エラー: --scan-kp だけでは実行できません（isBeam を既定で読まないため）。\n"
            "  次のどちらかを付けてください:\n"
            "    --close-dist-max MM   例: --close-dist-max 20\n"
            "    --read-isbeam         isBeam 枝が安定しているときのみ",
            file=sys.stderr,
        )
        sys.exit(3)
    be = "pyroot" if args.pyroot else args.backend
    set_path_helix_3d(args.rootfile, read_isbeam=args.read_isbeam, backend=be)
    if not _io_ready():
        sys.exit(1)
    entry = args.entry
    if args.scan_kp:
        found = scan_kp_candidate(args.start, args.max_scan, args.close_dist_max)
        if found is None:
            sys.exit(2)
        entry = found
    elif args.run is not None:
        if args.entry >= 0:
            print(
                f"注意: 位置引数 entry={args.entry} は --run 指定時は無視し、"
                f"run={args.run} を検索します。",
                flush=True,
            )
        fe = find_entry_helix_3d(
            args.run,
            args.event,
            start_entry=0,
            max_scan=args.max_run_scan,
        )
        if fe is None:
            sys.exit(4)
        entry = fe
    elif args.event is not None:
        if args.entry >= 0:
            print(
                f"注意: 位置引数 entry={args.entry} は --event のみ指定時は無視し、"
                f"event_number={args.event} を検索します。",
                flush=True,
            )
        fe = find_entry_by_event_number(
            args.event,
            start_entry=0,
            max_scan=args.max_run_scan,
        )
        if fe is None:
            sys.exit(4)
        entry = fe
    if args.browse and args.save is None:
        event_helix_3d_browse(first_entry=entry, kp_only=args.kp)
    else:
        event_helix_3d(
            entry,
            kp_only=args.kp,
            save_png=args.save,
            save_dpi=args.save_dpi,
        )


if __name__ == "__main__":
    main(sys.argv[1:])
