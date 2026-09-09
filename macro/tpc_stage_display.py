#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TPC 段階別 3D イベントディスプレイ（スライド用）。

DstTPCHelixTracking 出力（"tpc" ツリー）だけを使い、解析の段階
    Cluster (track_cluster_*_center) -> Helix (helix_*)
    -> PID (pid)  -> Vertex (vtxTpc, lambda/k0 mass)
を切り替えて描画する。座標変換・PID 色分け・検出器フレーム描画・軸ガイド等は
check_tpc_helix_track_3d.py（同じ macro/ 内）の実装をそのまま import して再利用する。
I/O は check_tpc_helix_track_3d.py と同じく既定で uproot（ROOT は使わない）。

段階の意味（各段階は前段までの内容を含む＝累積表示）:
  cluster : 各トラックに割り当てられたクラスタ重心の点群のみ（ニュートラル色）
  helix   : cluster に、ヘリックスフィット曲線を重ねる（まだニュートラル色）
  pid     : 点群・曲線の色を PID 色（K=blue, p=red, pi=green）に変更、凡例を表示
  vertex  : 頂点マーカー、closeDist / Lambda・K0 質量の HUD テキストを追加

現時点ではまだアニメーション（GIF/MP4）化はしていない。まずは各段階の静止画を
「いい感じに」表示できることを確認するためのスクリプト。

Usage（python から、check_tpc_helix_track_3d.py と同じ使用感）:
  python3
  >>> import tpc_stage_display as sd
  >>> sd.set_path("run02601_TPCHelix.root")
  >>> sd.show(-1, stage="helix")        # ランダムイベント、helix 段階まで表示
  >>> sd.show(12345, stage="vertex")    # entry 指定、全段階表示

CLI:
  python3 tpc_stage_display.py file.root --entry 12345 --stage pid --save out.png
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import warnings
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_tpc_helix_track_3d as base  # noqa: E402  座標変換・PID色・フレーム描画を再利用

STAGES = ["cluster", "helix", "pid", "vertex"]
# 初期カメラ視点。check_tpc_helix_track_3d.py の HELIX_VIEW_ELEV/AZIM (50.0/-90.0) は
# 上から見下ろす角度が強いため、本スクリプトではもう少し低い（水平に近い）角度にする。
VIEW_ELEV = 18.0
VIEW_AZIM = -90.0
NEUTRAL_COLOR = "0.35"  # PID 確定前のグレー（未使用、後方互換のため残置）
# PID 確定前（cluster/helix 段階）のトラック色。カラフルにせず単色（オレンジ）に統一する。
NEUTRAL_TRACK_COLORS = ["tab:orange"]
# helix 段階（stage="helix" のみ）でのヘリックス線の色。cluster 点（オレンジのまま）と区別しやすくする。
# pid/vertex 段階では従来通り PID 色を点・線の両方に使う（ここでは変更しない）。
HELIX_STAGE_LINE_COLOR = "black"
VERTEX_CLOSE_DIST_MAX_DEFAULT = 50.0  # mm; これより遠いペアの最近接点は表示しない（見た目のノイズ抑制）

# --- TPC pad geometry: tpc::padParameter（include/TPCPadHelper.hh）の写し ---
# columns: layerID, numOfPad, radius[mm], numOfDivision, dummy, length[mm]
PAD_PARAMETER = [
    (0, 48, 14.75, 48, 0, 9.0),
    (1, 48, 24.25, 48, 0, 9.0),
    (2, 72, 33.75, 72, 0, 9.0),
    (3, 96, 43.25, 96, 0, 9.0),
    (4, 120, 52.75, 120, 0, 9.0),
    (5, 144, 62.25, 144, 0, 9.0),
    (6, 168, 71.75, 168, 0, 9.0),
    (7, 192, 81.25, 192, 0, 9.0),
    (8, 216, 90.75, 216, 0, 9.0),
    (9, 240, 100.25, 240, 0, 9.0),
    (10, 208, 111.5, 241, 0, 12.5),
    (11, 218, 124.5, 271, 0, 12.5),
    (12, 230, 137.5, 300, 0, 12.5),
    (13, 214, 150.5, 330, 0, 12.5),
    (14, 212, 163.5, 360, 0, 12.5),
    (15, 214, 176.5, 390, 0, 12.5),
    (16, 220, 189.5, 420, 0, 12.5),
    (17, 224, 202.5, 449, 0, 12.5),
    (18, 232, 215.5, 479, 0, 12.5),
    (19, 238, 228.5, 509, 0, 12.5),
    (20, 244, 241.5, 539, 0, 12.5),
    (21, 232, 254.5, 569, 0, 12.5),
    (22, 218, 267.5, 599, 0, 12.5),
    (23, 210, 280.5, 628, 0, 12.5),
    (24, 206, 293.5, 658, 0, 12.5),
    (25, 202, 306.5, 688, 0, 12.5),
    (26, 200, 319.5, 718, 0, 12.5),
    (27, 196, 332.5, 748, 0, 12.5),
    (28, 178, 345.5, 777, 0, 12.5),
    (29, 130, 358.5, 807, 0, 12.5),
    (30, 108, 371.5, 837, 0, 12.5),
    (31, 90, 384.5, 867, 0, 12.5),
]
NUM_LAYERS_TPC = len(PAD_PARAMETER)  # 32, tpc::NumOfLayersTPC
# パッド平面を表示するドリフト方向(局所 y)の位置。実際のヒット y はドリフト時間から別途決まるため、
# ここは「読み出し面がある」ことを示す視覚表現としてフレーム下端に固定するだけ（物理的厳密位置ではない）。
PAD_PLANE_Y = -310.0


def _pad_center_theta_deg(layer: int, row: float) -> float:
    """tpc::GetTheta(layer, m_row) と同じ（パッド中心の方位角 [deg]）。"""
    _, n_pad, _, n_div, _, _ = PAD_PARAMETER[layer]
    d_theta = 360.0 / n_div
    s_theta = 180.0 - d_theta * n_pad / 2.0
    return s_theta + (row + 0.5) * d_theta - 180.0


def _pad_boundary_xz(layer: int, row: int):
    """パッド (layer, row) の 4 頂点（局所 x, z）。中心角 ± 半分の扇形近似（tpc::InitializeHistograms と同じ形）。"""
    _, _, radius, n_div, _, length = PAD_PARAMETER[layer]
    d_theta = 360.0 / n_div
    theta_c = _pad_center_theta_deg(layer, row)
    theta1 = math.radians(theta_c - d_theta / 2.0)
    theta2 = math.radians(theta_c + d_theta / 2.0)
    r_min = radius - length / 2.0
    r_max = radius + length / 2.0
    xs = (r_max * math.sin(theta1), r_max * math.sin(theta2), r_min * math.sin(theta2), r_min * math.sin(theta1))
    zs = (
        r_max * math.cos(theta1) + base.TPC_Z_TARGET,
        r_max * math.cos(theta2) + base.TPC_Z_TARGET,
        r_min * math.cos(theta2) + base.TPC_Z_TARGET,
        r_min * math.cos(theta1) + base.TPC_Z_TARGET,
    )
    return xs, zs


def _build_all_pad_polygons():
    """
    全パッド（32層、計 5768 枚）の表示座標系ポリゴン頂点リストと、(layer,row) -> polys 内 index の対応を作る。
    モジュール import 時に一度だけ計算してキャッシュする（_TPC_FRAME_SEGMENTS と同じパターン）。
    """
    polys: List[List[List[float]]] = []
    index: dict = {}
    for layer in range(NUM_LAYERS_TPC):
        n_pad = int(PAD_PARAMETER[layer][1])
        for row in range(n_pad):
            xs, zs = _pad_boundary_xz(layer, row)
            poly = [list(base.tpc_local_to_display(x, PAD_PLANE_Y, z)) for x, z in zip(xs, zs)]
            index[(layer, row)] = len(polys)
            polys.append(poly)
    return polys, index


_PAD_POLY_VERTS, _PAD_POLY_INDEX = _build_all_pad_polygons()


@dataclass
class StageExtra:
    """base.EventDataHelix3D には無い、段階表示専用の追加ブランチ。"""

    has_cluster: bool = False
    track_cluster_x: List[List[float]] = field(default_factory=list)
    track_cluster_y: List[List[float]] = field(default_factory=list)
    track_cluster_z: List[List[float]] = field(default_factory=list)
    has_hitpad: bool = False
    hitlayer: List[List[float]] = field(default_factory=list)
    track_cluster_row_center: List[List[float]] = field(default_factory=list)
    has_lambda: bool = False
    lambda_mass: List[float] = field(default_factory=list)
    lambda_vtx_x: List[float] = field(default_factory=list)
    lambda_vtx_y: List[float] = field(default_factory=list)
    lambda_vtx_z: List[float] = field(default_factory=list)
    has_k0: bool = False
    k0_mass: List[float] = field(default_factory=list)
    k0_vtx_x: List[float] = field(default_factory=list)
    k0_vtx_y: List[float] = field(default_factory=list)
    k0_vtx_z: List[float] = field(default_factory=list)


g_extra = StageExtra()
g_fig = None
g_ax = None
g_checkbox_ax = None
g_checkbox = None
g_track_artists: dict = {}  # itrack -> list[Artist]（interactive 表示切替用）
g_track_visible: dict = {}  # itrack -> bool（interactive でのチェックボックス状態）
g_vertex_artists: List = []  # vertex 段階の Artist 一覧（インデックス = interactive 表示切替用の vtx{i}）
g_vertex_visible: dict = {}  # ivtx -> bool（interactive でのチェックボックス状態）
g_current_render_kwargs: dict = {}  # export_current() が再描画に使う render_stage 引数


def set_path(path: str, backend: str = "uproot") -> None:
    """base 側の set_path_helix_3d をそのまま使う（既定 uproot）。"""
    base.set_path_helix_3d(path, backend=backend)


def _load_extra_uproot(entry: int) -> StageExtra:
    t = base.g_uproot_tree
    ex = StageExtra()
    ex.has_cluster = base._uproot_has_branch(t, "track_cluster_x_center")
    ex.has_hitpad = base._uproot_has_branch(t, "hitlayer") and base._uproot_has_branch(t, "track_cluster_row_center")
    ex.has_lambda = base._uproot_has_branch(t, "lambda_mass")
    ex.has_k0 = base._uproot_has_branch(t, "k0_mass")

    names: List[str] = []
    if ex.has_cluster:
        names += ["track_cluster_x_center", "track_cluster_y_center", "track_cluster_z_center"]
    if ex.has_hitpad:
        names += ["hitlayer", "track_cluster_row_center"]
    if ex.has_lambda:
        names += ["lambda_mass", "lambda_vtx_x", "lambda_vtx_y", "lambda_vtx_z"]
    if ex.has_k0:
        names += ["k0_mass", "k0_vtx_x", "k0_vtx_y", "k0_vtx_z"]
    if not names:
        return ex

    chunk = t.arrays(names, library="ak", entry_start=entry, entry_stop=entry + 1)
    nt = base.g_event_helix.ntTpc
    if ex.has_cluster:
        ex.track_cluster_x = base._cell0_nested_float(chunk, "track_cluster_x_center", nt)
        ex.track_cluster_y = base._cell0_nested_float(chunk, "track_cluster_y_center", nt)
        ex.track_cluster_z = base._cell0_nested_float(chunk, "track_cluster_z_center", nt)
    if ex.has_hitpad:
        ex.hitlayer = base._cell0_nested_float(chunk, "hitlayer", nt)
        ex.track_cluster_row_center = base._cell0_nested_float(chunk, "track_cluster_row_center", nt)
    if ex.has_lambda:
        ex.lambda_mass = base._cell0_vec_float(chunk, "lambda_mass")
        ex.lambda_vtx_x = base._cell0_vec_float(chunk, "lambda_vtx_x")
        ex.lambda_vtx_y = base._cell0_vec_float(chunk, "lambda_vtx_y")
        ex.lambda_vtx_z = base._cell0_vec_float(chunk, "lambda_vtx_z")
    if ex.has_k0:
        ex.k0_mass = base._cell0_vec_float(chunk, "k0_mass")
        ex.k0_vtx_x = base._cell0_vec_float(chunk, "k0_vtx_x")
        ex.k0_vtx_y = base._cell0_vec_float(chunk, "k0_vtx_y")
        ex.k0_vtx_z = base._cell0_vec_float(chunk, "k0_vtx_z")
    return ex


def _load_extra_pyroot(entry: int) -> StageExtra:
    t = base.g_tree_helix
    ex = StageExtra()
    ex.has_cluster = bool(t.GetBranch("track_cluster_x_center"))
    ex.has_hitpad = bool(t.GetBranch("hitlayer")) and bool(t.GetBranch("track_cluster_row_center"))
    ex.has_lambda = bool(t.GetBranch("lambda_mass"))
    ex.has_k0 = bool(t.GetBranch("k0_mass"))

    nt = base.g_event_helix.ntTpc
    if ex.has_cluster:
        ex.track_cluster_x = base._nested_vec_double_from_tree(t, "track_cluster_x_center", nt)
        ex.track_cluster_y = base._nested_vec_double_from_tree(t, "track_cluster_y_center", nt)
        ex.track_cluster_z = base._nested_vec_double_from_tree(t, "track_cluster_z_center", nt)
    if ex.has_hitpad:
        ex.hitlayer = base._nested_vec_double_from_tree(t, "hitlayer", nt)
        ex.track_cluster_row_center = base._nested_vec_double_from_tree(t, "track_cluster_row_center", nt)
    if ex.has_lambda:
        ex.lambda_mass = base._vec_double_to_list(t.lambda_mass)
        ex.lambda_vtx_x = base._vec_double_to_list(t.lambda_vtx_x)
        ex.lambda_vtx_y = base._vec_double_to_list(t.lambda_vtx_y)
        ex.lambda_vtx_z = base._vec_double_to_list(t.lambda_vtx_z)
    if ex.has_k0:
        ex.k0_mass = base._vec_double_to_list(t.k0_mass)
        ex.k0_vtx_x = base._vec_double_to_list(t.k0_vtx_x)
        ex.k0_vtx_y = base._vec_double_to_list(t.k0_vtx_y)
        ex.k0_vtx_z = base._vec_double_to_list(t.k0_vtx_z)
    return ex


def load_event(entry: int = -1) -> Optional[int]:
    """base.load_event_helix_3d に加え、cluster/lambda/k0 の追加ブランチを読む。"""
    global g_extra
    res = base.load_event_helix_3d(entry)
    if res is None:
        return None
    g_extra = _load_extra_uproot(res) if base.g_backend == "uproot" else _load_extra_pyroot(res)
    return res


def _symmetrize_about_origin(xmin, xmax, ymin, ymax, zmin, zmax):
    """
    表示座標系の原点（= TPC 機械中心, フレームの幾何中心）を軸範囲の中心に揃える。
    matplotlib の 3D 回転はマウスドラッグ操作の中心ではなく xlim/ylim/zlim の中心を軸に
    行われるため、ヒット分布に偏って非対称な範囲だと「TPC の真ん中で回っていない」ように見える。
    """
    def sym(lo, hi):
        r = max(abs(lo), abs(hi), 1e-6)
        return -r, r

    return (*sym(xmin, xmax), *sym(ymin, ymax), *sym(zmin, zmax))


def _collect_bounds_with_clusters(ev, ex: StageExtra):
    """base._collect_bounds（hit/helix ベース）に、cluster 点群の範囲も加える。"""
    xmin, xmax, ymin, ymax, zmin, zmax = base._collect_bounds(ev)
    if not ex.has_cluster:
        return xmin, xmax, ymin, ymax, zmin, zmax

    xs, ys, zs = [], [], []
    for i in range(ev.ntTpc):
        if i >= len(ex.track_cluster_x):
            continue
        for j in range(len(ex.track_cluster_x[i])):
            mx, my, mz = base.tpc_local_to_display(
                ex.track_cluster_x[i][j], ex.track_cluster_y[i][j], ex.track_cluster_z[i][j]
            )
            xs.append(mx)
            ys.append(my)
            zs.append(mz)
    if xs:
        xmin, xmax = min(xmin, min(xs)), max(xmax, max(xs))
        ymin, ymax = min(ymin, min(ys)), max(ymax, max(ys))
        zmin, zmax = min(zmin, min(zs)), max(zmax, max(zs))
    return xmin, xmax, ymin, ymax, zmin, zmax


def _draw_vertices(ax, ev, ex: StageExtra, close_dist_max: float) -> List:
    """
    頂点マーカーをまとめて描画する。Lambda/K0 等の粒子種は区別せず、すべて同じ
    見た目（マゼンタの X）の汎用 "vertex" として描画する。
    凡例には、分かる範囲でどのトラックの組み合わせかを付記する:
      - 全トラックペアの最近接点（vtxTpc[i][j], closeDistTpc[i][j]）のうち closeDist <= close_dist_max のもの
        （C++/旧 Python 版の "ntTpc==2 のみ" という制約はここでは外している。vtxTpc は元々
        [it][it_pair] の全ペア行列: DstTPCHelixTracking.cc FillHelixPairKinematics）。
        トラック番号 i, j が既知なので "vertex (tr{i}-tr{j})" とする。
      - Lambda / K0 の崩壊点（lambda_vtx_*, k0_vtx_*）。現状の DST 出力にはどのトラックの組か
        を示すブランチが無い（track_id 系ブランチ非搭載）ため、組み合わせは表示できず、
        まとめて先頭の 1 点にのみ "vertex" と付記する（凡例の重複を避けるため）。
    戻り値: 描画した Artist のリスト（インデックスが interactive チェックボックスの vtx{i} に対応）。
    """
    artists: List = []
    generic_labeled = False

    def _scatter(vx: float, vy: float, vz: float, label: Optional[str]) -> None:
        nonlocal generic_labeled
        if not (math.isfinite(vx) and math.isfinite(vy) and math.isfinite(vz)):
            return
        mx, my, mz = base.tpc_local_to_display(vx, vy, vz)
        if label is not None:
            lbl = label
        elif not generic_labeled:
            lbl = "vertex"
            generic_labeled = True
        else:
            lbl = None
        sc = ax.scatter([mx], [my], [mz], c="magenta", marker="x", s=140, linewidths=2.2, label=lbl, zorder=9)
        artists.append(sc)

    if len(ev.vtxTpc) != 0:
        for i in range(ev.ntTpc):
            if i >= len(ev.vtxTpc):
                continue
            for j in range(i + 1, ev.ntTpc):
                if j >= len(ev.vtxTpc[i]):
                    continue
                cd = ev.closeDistTpc[i][j] if (i < len(ev.closeDistTpc) and j < len(ev.closeDistTpc[i])) else float("nan")
                if math.isfinite(cd) and cd > close_dist_max:
                    continue
                _scatter(ev.vtxTpc[i][j], ev.vtyTpc[i][j], ev.vtzTpc[i][j], label=f"vertex (tr{i}-tr{j})")

    for k in range(len(ex.lambda_vtx_x)):
        _scatter(ex.lambda_vtx_x[k], ex.lambda_vtx_y[k], ex.lambda_vtx_z[k], label=None)
    for k in range(len(ex.k0_vtx_x)):
        _scatter(ex.k0_vtx_x[k], ex.k0_vtx_y[k], ex.k0_vtx_z[k], label=None)

    return artists


def _draw_pad_background(ax) -> None:
    """全パッド（32層、計 5768 枚）を薄いグレーの扇形として描画する（draw_pads=True 時の背景）。"""
    coll = Poly3DCollection(
        _PAD_POLY_VERTS, facecolor="0.88", edgecolor="0.75", linewidths=0.15, alpha=0.35, zorder=1,
    )
    ax.add_collection3d(coll)


def _draw_hit_pads_for_track(ax, ex: StageExtra, itrack: int, color, visible: bool):
    """
    hitlayer[itrack] / track_cluster_row_center[itrack]（同じインデックスで対応:
    DstTPCHelixTracking.cc の ProcessOneTrackHit）から、そのトラックがヒットしたパッドを
    color で塗った 1 つの Poly3DCollection として返す（ヒットが無ければ None）。
    """
    if not ex.has_hitpad or itrack >= len(ex.hitlayer) or itrack >= len(ex.track_cluster_row_center):
        return None
    layers = ex.hitlayer[itrack]
    rows = ex.track_cluster_row_center[itrack]
    polys = []
    for lay, rw in zip(layers, rows):
        if not (math.isfinite(lay) and math.isfinite(rw)):
            continue
        key = (int(round(lay)), int(round(rw)))
        idx = _PAD_POLY_INDEX.get(key)
        if idx is not None:
            polys.append(_PAD_POLY_VERTS[idx])
    if not polys:
        return None
    coll = Poly3DCollection(
        polys, facecolor=color, edgecolor=color, linewidths=1.4, alpha=1.0, zorder=3, visible=visible,
    )
    ax.add_collection3d(coll)
    return coll


def render_stage(
    ax,
    ev,
    ex: StageExtra,
    stage: str = "vertex",
    entry_label=None,
    vertex_close_dist_max: float = VERTEX_CLOSE_DIST_MAX_DEFAULT,
    track_visible: Optional[dict] = None,
    vertex_visible: Optional[dict] = None,
    draw_pads: bool = False,
    for_export: bool = False,
):
    """
    指定 stage までの要素を累積描画する（内部で ax.clear() する）。
    track_visible: {itrack: bool} で特定トラックを非表示にできる（省略時は全トラック表示）。
    vertex_visible: {ivtx: bool} で特定の頂点マーカーを個別に非表示にできる（省略時は全頂点表示）。
      ivtx は _draw_vertices() が描画する順（トラックペア最近接点 -> Lambda -> K0）のインデックス。
    draw_pads: True で全 TPC パッド（背景・薄灰）とヒットパッド（トラック色）を描画する。
      パッド数が多く（32層・計5768枚）描画が重くなるため既定は False。
    for_export: True なら凡例（トラック名・頂点ラベル等）を描画しない
      （画像保存時に画面を占有する凡例を省くため。show(save=...) から自動的に True になる）。
    戻り値: ({itrack: [Artist, ...]}, [Artist, ...]) のタプル
      （それぞれ interactive なトラック / 頂点ごとの表示切替に使う。後者のリストの
      インデックスが vertex_visible のキーおよび interactive チェックボックスの vtx{i} に対応）。
    """
    if stage not in STAGES:
        raise ValueError(f"unknown stage {stage!r}; must be one of {STAGES}")
    stage_idx = STAGES.index(stage)
    show_helix = stage_idx >= STAGES.index("helix")
    show_pid = stage_idx >= STAGES.index("pid")
    show_vertex = stage_idx >= STAGES.index("vertex")
    track_visible = track_visible or {}
    vertex_visible = vertex_visible or {}

    ax.clear()
    xmin, xmax, ymin, ymax, zmin, zmax = _collect_bounds_with_clusters(ev, ex)
    xmin, xmax, ymin, ymax, zmin, zmax = base._merge_span_with_frame(xmin, xmax, ymin, ymax, zmin, zmax)
    xmin, xmax, ymin, ymax, zmin, zmax = _symmetrize_about_origin(xmin, xmax, ymin, ymax, zmin, zmax)

    base._draw_tpc_frame(ax)
    if draw_pads:
        _draw_pad_background(ax)

    track_artists: dict = {}
    for itrack in range(ev.ntTpc):
        pid_code = ev.pid[itrack] if itrack < len(ev.pid) else 0
        if show_pid:
            color, pk = base.pid_display_color(pid_code)
            lbl = f"tr{itrack} {pk} (pid=0x{pid_code:x}, {base.decode_pid_candidates(pid_code)})"
        else:
            color = NEUTRAL_TRACK_COLORS[itrack % len(NEUTRAL_TRACK_COLORS)]
            lbl = None

        artists: List = []
        visible = track_visible.get(itrack, True)

        if ex.has_cluster and itrack < len(ex.track_cluster_x) and ex.track_cluster_x[itrack]:
            cx = np.asarray(ex.track_cluster_x[itrack], dtype=float)
            cy = np.asarray(ex.track_cluster_y[itrack], dtype=float)
            cz = np.asarray(ex.track_cluster_z[itrack], dtype=float)
            mx, my, mz = base.tpc_local_to_display_vec(cx, cy, cz)
            # helix 曲線が乗る段階（helix 以降）はヘリックス線を主役にし、点は控えめにする。
            cluster_marker_size = 8 if show_helix else 16
            sc = ax.scatter(mx, my, mz, color=color, marker="o", s=cluster_marker_size, alpha=0.85, visible=visible)
            artists.append(sc)

        if show_helix:
            cxh = ev.helix_cx[itrack] if itrack < len(ev.helix_cx) else 0.0
            cyh = ev.helix_cy[itrack] if itrack < len(ev.helix_cy) else 0.0
            z0 = ev.helix_z0[itrack] if itrack < len(ev.helix_z0) else 0.0
            rr = ev.helix_r[itrack] if itrack < len(ev.helix_r) else 0.0
            dz = ev.helix_dz[itrack] if itrack < len(ev.helix_dz) else 0.0
            tr = base.helix_draw_theta_range(ev, itrack)
            if tr is not None and math.isfinite(float(rr)):
                t0, t1 = tr
                nt = max(8, base.helix_polyline_sample_count(t0, t1))
                thetas = np.linspace(t0, t1, nt)
                xh, yh, zh = base.helix_xyz(cxh, cyh, z0, rr, dz, thetas)
                mx, my, mz = base.tpc_local_to_display_vec(xh, yh, zh)
                fin = np.isfinite(mx) & np.isfinite(my) & np.isfinite(mz)
                if np.any(fin):
                    line_color = HELIX_STAGE_LINE_COLOR if stage == "helix" else color
                    (ln,) = ax.plot(
                        mx[fin], my[fin], mz[fin],
                        color=line_color, linewidth=2.2, linestyle=":", label=lbl, zorder=6,
                        visible=visible,
                    )
                    artists.append(ln)

        if draw_pads:
            pad_coll = _draw_hit_pads_for_track(ax, ex, itrack, color, visible)
            if pad_coll is not None:
                artists.append(pad_coll)

        track_artists[itrack] = artists

    vertex_artists: List = []
    if show_vertex:
        vertex_artists = _draw_vertices(ax, ev, ex, vertex_close_dist_max)
        for i, artist in enumerate(vertex_artists):
            artist.set_visible(vertex_visible.get(i, True))

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)

    if not for_export:
        title = f"TPC [{stage}] (Run {ev.runnum}, Event {ev.evnum}"
        if entry_label is not None:
            title += f", entry {entry_label}"
        title += ")"
        ax.set_title(title, fontsize=10)
    if show_pid and not for_export:
        ax.legend(loc="upper left", fontsize=7)

    if show_vertex:
        hud_lines = []
        if ex.has_lambda and ex.lambda_mass:
            vals = ", ".join(f"{m:.4g}" for m in ex.lambda_mass if math.isfinite(m))
            if vals:
                hud_lines.append(f"Lambda mass: {vals} GeV/c^2")
        if ex.has_k0 and ex.k0_mass:
            vals = ", ".join(f"{m:.4g}" for m in ex.k0_mass if math.isfinite(m))
            if vals:
                hud_lines.append(f"K0 mass: {vals} GeV/c^2")
        if hud_lines:
            ax.text2D(
                0.02, 0.02, "\n".join(hud_lines), transform=ax.transAxes, fontsize=8,
                va="bottom", ha="left",
                bbox=dict(boxstyle="round", fc="white", ec="0.5", alpha=0.85),
            )

    base._style_ax_without_builtin_axes(ax)
    base._draw_manual_xyz_axes(ax, xmin, xmax, ymin, ymax, zmin, zmax)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    ax.view_init(elev=VIEW_ELEV, azim=VIEW_AZIM)
    with warnings.catch_warnings():
        # interactive 版は checkbox 用 Axes を add_axes で手動配置しているため
        # tight_layout 非対応の UserWarning が出るが、実害はない。
        warnings.simplefilter("ignore", UserWarning)
        plt.tight_layout()
    return track_artists, vertex_artists


g_interactive_mode = False


def _ensure_figure(interactive: bool = False) -> None:
    """
    interactive=True なら、3D Axes の右にトラック ON/OFF チェックボックス用の Axes を確保する。
    既存 figure の interactive 状態が要求と異なる場合は作り直す。
    """
    global g_fig, g_ax, g_checkbox_ax, g_interactive_mode
    if g_fig is not None and plt.fignum_exists(g_fig.number) and g_interactive_mode == interactive:
        return
    if g_fig is not None and plt.fignum_exists(g_fig.number):
        plt.close(g_fig)
    g_interactive_mode = interactive
    if interactive:
        g_fig = plt.figure(figsize=(10, 8))
        g_ax = g_fig.add_axes([0.03, 0.03, 0.74, 0.94], projection="3d")
        g_checkbox_ax = g_fig.add_axes([0.80, 0.30, 0.18, 0.4])
    else:
        g_fig = plt.figure(figsize=(8, 8))
        g_ax = g_fig.add_subplot(111, projection="3d")
        g_checkbox_ax = None


def _setup_track_checkboxes(ev, track_artists: dict, vertex_artists: List) -> None:
    """トラックごと + vertex 全体の ON/OFF チェックボックスを g_checkbox_ax に配置する（interactive 表示専用）。"""
    global g_checkbox, g_track_artists, g_track_visible, g_vertex_artists, g_vertex_visible
    g_track_artists = track_artists
    g_track_visible = {itrack: True for itrack in range(ev.ntTpc)}
    g_vertex_artists = vertex_artists
    g_vertex_visible = {i: True for i in range(len(vertex_artists))}
    if g_checkbox_ax is None:
        return
    g_checkbox_ax.clear()
    g_checkbox_ax.set_title("tracks/vertices\n('e': export PNG)", fontsize=9)
    labels = [f"tr{itrack}" for itrack in range(ev.ntTpc)]
    labels += [f"vtx{i}" for i in range(len(vertex_artists))]
    if not labels:
        g_checkbox = None
        return
    g_checkbox = CheckButtons(g_checkbox_ax, labels, [True] * len(labels))
    g_checkbox.on_clicked(_on_checkbox_clicked)


def _on_checkbox_clicked(label: str) -> None:
    if label.startswith("vtx"):
        ivtx = int(label[3:])  # "vtx2" -> 2
        g_vertex_visible[ivtx] = not g_vertex_visible.get(ivtx, True)
        if 0 <= ivtx < len(g_vertex_artists):
            g_vertex_artists[ivtx].set_visible(g_vertex_visible[ivtx])
    else:
        itrack = int(label[2:])  # "tr3" -> 3
        g_track_visible[itrack] = not g_track_visible.get(itrack, True)
        for artist in g_track_artists.get(itrack, []):
            artist.set_visible(g_track_visible[itrack])
    if g_fig is not None:
        g_fig.canvas.draw_idle()


def export_current(path: Optional[str] = None, dpi: int = 300) -> Optional[str]:
    """
    直近の show(interactive=True) で表示中のイベントを、その時点でチェックボックスで
    非表示にしたトラック・頂点の状態を保ったまま、タイトル・凡例なし／背景透過のエクスポート用
    画像として別ファイルに保存する（画面のウィンドウ自体はそのまま維持される）。
    path 省略時は "tpc_export_run{run}_ev{event}_{stage}.png" を使う。
    """
    if not g_current_render_kwargs:
        print("Error: 表示中のイベントがありません。先に show(..., interactive=True) を呼んでください。")
        return None
    ev = base.g_event_helix
    if path is None:
        path = f"tpc_export_run{ev.runnum}_ev{ev.evnum}_{g_current_render_kwargs['stage']}.png"
    export_fig = plt.figure(figsize=(8, 8))
    export_ax = export_fig.add_subplot(111, projection="3d")
    render_stage(
        export_ax, ev, g_extra,
        track_visible=dict(g_track_visible),
        vertex_visible=dict(g_vertex_visible),
        for_export=True,
        **g_current_render_kwargs,
    )
    base._helix_apply_transparent_bg(export_fig, export_ax)
    export_fig.savefig(path, dpi=dpi, transparent=True)
    plt.close(export_fig)
    print(f"Saved: {path}")
    return path


def _on_interactive_key_press(event) -> None:
    # 's' は matplotlib 既定のショートカット（Figure を名前を付けて保存）と衝突して
    # export_current() が発火しないことがあるため 'e'（export）を使う。
    if getattr(event, "key", "") not in ("e", "E"):
        return
    export_current()


def show(
    entry: int = -1,
    stage: str = "vertex",
    save: Optional[str] = None,
    save_dpi: int = 300,
    interactive: bool = False,
    vertex_close_dist_max: float = VERTEX_CLOSE_DIST_MAX_DEFAULT,
    draw_pads: bool = False,
) -> Optional[int]:
    """
    1 イベントを読み込み、指定 stage まで累積描画する。
    save 指定時は PNG 保存のみ（画面表示なし、interactive は無視される）。
    保存画像はスライド等への貼り付け用に、タイトル・凡例なし、背景透過、高解像度（既定 dpi=300）にする。
    interactive=True なら、トラックごと + vertex（stage="vertex" の場合のみ）の
    ON/OFF チェックボックスを画面に表示する
    （X11 転送など、実際に matplotlib ウィンドウが表示できる環境が必要）。
    チェックボックスで好きなトラック/vertex を非表示にしたあと、その状態を PNG に書き出すには
    二通りある:
      (a) 表示中のウィンドウで 'e' キーを押す（'s' は matplotlib 標準の保存ダイアログと
          衝突するため使わない）。
      (b) ウィンドウを閉じたあと、Python 側で export_current("out.png") を呼ぶ
          （チェックボックスの状態は show() を呼び直すまで記憶されている）。
    (a) が反応しない環境では (b) を使ってください。
    draw_pads=True で全 TPC パッド＋ヒットパッドを描画する（重いので既定 False）。
    """
    global g_current_render_kwargs
    res = load_event(entry)
    if res is None:
        return None
    use_interactive = bool(interactive) and not save
    _ensure_figure(interactive=use_interactive)
    g_current_render_kwargs = dict(
        stage=stage, entry_label=res,
        vertex_close_dist_max=vertex_close_dist_max, draw_pads=draw_pads,
    )
    track_artists, vertex_artists = render_stage(
        g_ax, base.g_event_helix, g_extra, for_export=bool(save), **g_current_render_kwargs,
    )
    if use_interactive:
        _setup_track_checkboxes(base.g_event_helix, track_artists, vertex_artists)
        g_fig.canvas.mpl_connect("key_press_event", _on_interactive_key_press)
    if save:
        base._helix_apply_transparent_bg(g_fig, g_ax)
        g_fig.savefig(save, dpi=save_dpi, transparent=True)
        print(f"Saved: {save}")
    else:
        if use_interactive:
            print("Tip: チェックボックスでトラック/vertex を非表示にした後、"
                  "'e' キーで現在の状態を PNG (透過・凡例なし) にエクスポートできます。")
            print("     反応しない場合は、ウィンドウを閉じてから "
                  "export_current(\"out.png\") を呼んでください（状態は記憶されています）。")
        plt.show()
    return res


def _parse_cli(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="TPC 段階別 3D イベントディスプレイ（Cluster/Helix/PID/Vertex）")
    p.add_argument("rootfile", help="DstTPCHelixTracking 出力 ROOT ファイル（tree 'tpc'）")
    p.add_argument("entry", type=int, nargs="?", default=-1, help="entry 番号（省略時 or 負値でランダム）")
    p.add_argument("--stage", choices=STAGES, default="vertex", help="表示する段階（既定 vertex=全段階）")
    p.add_argument("--save", help="PNG 保存先パス（省略時は画面表示）")
    p.add_argument("--save-dpi", type=int, default=300)
    p.add_argument("--backend", choices=["auto", "uproot", "pyroot"], default="uproot")
    p.add_argument(
        "--interactive", action="store_true",
        help="トラック ON/OFF チェックボックス付きで画面表示（--save 指定時は無視）。X11 転送等が必要。",
    )
    p.add_argument(
        "--vertex-close-dist-max", type=float, default=VERTEX_CLOSE_DIST_MAX_DEFAULT,
        help=f"vertex 段階でトラックペアの最近接点を描く closeDist 上限 [mm]（既定 {VERTEX_CLOSE_DIST_MAX_DEFAULT}）",
    )
    p.add_argument(
        "--draw-pads", action="store_true",
        help="全 TPC パッド（32層・計5768枚、薄灰）とヒットパッド（トラック色）を描画する（重いので既定オフ）。",
    )
    return p.parse_args(argv)


def main(argv: List[str]) -> None:
    args = _parse_cli(argv)
    set_path(args.rootfile, backend=args.backend)
    if not base._io_ready():
        sys.exit(1)
    res = show(
        args.entry, stage=args.stage, save=args.save, save_dpi=args.save_dpi,
        interactive=args.interactive, vertex_close_dist_max=args.vertex_close_dist_max,
        draw_pads=args.draw_pads,
    )
    if res is None:
        sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1:])
