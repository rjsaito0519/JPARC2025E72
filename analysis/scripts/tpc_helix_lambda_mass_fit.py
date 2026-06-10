#!/usr/bin/env python3
"""``lambda_mass`` ヒストグラムに Gauss + 線形背景をフィットし、シグナル yield を出す。

縦軸は **counts / bin**。連続モデルは **counts / (GeV/c^2)**（密度）とし、
フィット・描画では ``density * bin_width`` をビン内容と比較する。

  N_signal = area = ∫ (Gauss の密度) dx   [events]

``area`` をビン内容と素直に比べると、ピーク高さ ~800 に対して area ~10 のように
**counts×GeV/c^2** になってしまう（次元の取り違え）。

依存: pip install uproot awkward matplotlib lmfit
使い方: python3 tpc_helix_lambda_mass_fit.py path/to/file.root
"""

from __future__ import annotations

import os
import sys

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot
from lmfit import Model

# --- 必要なら編集（tpc_helix_lambda_mass.py と揃える） ---
TREE = "tpc"
BINS = 400
XMIN, XMAX = 1.05, 1.25  # GeV/c^2
FIT_XMIN, FIT_XMAX = 1.09, 1.12  # フィット範囲 [GeV/c^2]
# ---------------------------------------------------------


def gauss_density(x, area, center, sigma):
    """Gauss 密度 [counts / (GeV/c^2)]。area = ∫密度 dx [events]。"""
    amp = area / (sigma * np.sqrt(2.0 * np.pi))
    return amp * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def linear_density(x, center, slope, intercept):
    """線形背景密度 [counts / (GeV/c^2)]。intercept は x=center での密度。"""
    return intercept + slope * (x - center)


def gauss_plus_linear_density(x, area, center, sigma, slope, intercept):
    return gauss_density(x, area, center, sigma) + linear_density(
        x, center, slope, intercept
    )


def gauss_plus_linear_bin(x, bin_width, area, center, sigma, slope, intercept):
    """ビンあたりの期待値 [counts / bin]。"""
    return gauss_plus_linear_density(x, area, center, sigma, slope, intercept) * bin_width


def yield_stderr_from_cov(fit_result, param_name: str = "area") -> float:
    err = fit_result.params[param_name].stderr
    return float(err) if err is not None else float("nan")


def integrate_density(x: np.ndarray, y_density: np.ndarray) -> float:
    """∫ (counts/(GeV/c^2)) dx → events。"""
    if len(x) < 2:
        return 0.0
    return float(np.sum(0.5 * (y_density[1:] + y_density[:-1]) * np.diff(x)))


def integrate_linear_in_window(
    center: float, slope: float, intercept: float, xmin: float, xmax: float
) -> float:
    """線形背景密度を [xmin, xmax] で積分 → events。"""
    width = xmax - xmin
    moment = 0.5 * ((xmax - center) ** 2 - (xmin - center) ** 2)
    return intercept * width + slope * moment


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} path/to/file.root")

    rootfile = sys.argv[1]

    with uproot.open(rootfile) as src:
        if TREE not in src:
            sys.exit(f"error: tree {TREE!r} not in {rootfile!r}")
        tree = src[TREE]
        if "lambda_mass" not in tree:
            sys.exit(
                "error: branch 'lambda_mass' missing "
                "(EnableReconstructLambda ビルドの ROOT を使う)"
            )
        jagged = tree["lambda_mass"].array(library="ak")

    flat = ak.flatten(jagged, axis=None)
    if len(flat) == 0:
        sys.exit("error: no lambda_mass entries")

    masses = np.asarray(ak.to_numpy(flat), dtype=np.float64)

    counts, edges = np.histogram(masses, bins=BINS, range=(XMIN, XMAX))
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = float(edges[1] - edges[0])

    fit_mask = (centers >= FIT_XMIN) & (centers <= FIT_XMAX)
    x_fit = centers[fit_mask]
    y_fit = counts[fit_mask].astype(np.float64)
    yerr = np.sqrt(np.maximum(y_fit, 1.0))

    # フィット窓内のデータ（S+B）生カウント
    n_data_bins = float(np.sum(y_fit))
    n_data_entries = int(np.count_nonzero((masses >= FIT_XMIN) & (masses <= FIT_XMAX)))

    model = Model(gauss_plus_linear_bin)
    peak_i = int(np.argmax(y_fit))
    params = model.make_params(
        area=float(np.sum(y_fit) * 0.3),
        center=float(x_fit[peak_i]),
        sigma=0.002,
        slope=0.0,
        intercept=float(np.median(y_fit) / bin_width),
        bin_width=bin_width,
    )
    params["bin_width"].set(vary=False)
    params["center"].set(value=1.115, min=FIT_XMIN, max=FIT_XMAX)
    params["sigma"].set(min=bin_width, max=0.02)
    params["area"].set(min=0)

    result = model.fit(y_fit, params, x=x_fit, weights=1.0 / yerr)

    p = result.best_values
    yield_full = p["area"]
    yield_full_err = yield_stderr_from_cov(result, "area")

    x_num = np.linspace(FIT_XMIN, FIT_XMAX, 2000)
    g_den = gauss_density(x_num, p["area"], p["center"], p["sigma"])
    s_fit = integrate_density(x_num, g_den)
    b_fit = integrate_linear_in_window(
        p["center"], p["slope"], p["intercept"], FIT_XMIN, FIT_XMAX
    )
    s_fit_err = (
        yield_full_err * (s_fit / yield_full)
        if yield_full > 0 and np.isfinite(yield_full_err)
        else float("nan")
    )

    peak_density = yield_full / (p["sigma"] * np.sqrt(2.0 * np.pi))
    peak_counts_pred = peak_density * bin_width

    sn_ratio = s_fit / b_fit if b_fit > 0 else float("nan")
    significance_b = s_fit / np.sqrt(b_fit) if b_fit > 0 else float("nan")
    significance_tot = (
        s_fit / np.sqrt(n_data_bins) if n_data_bins > 0 else float("nan")
    )

    print(f"bin width: {bin_width:.6f} GeV/c^2")
    print(f"fit range: [{FIT_XMIN}, {FIT_XMAX}] GeV/c^2")
    print(result.fit_report())
    print()
    print("--- counts in fit window (data) ---")
    print(f"  N_data (sum of histogram bins):  {n_data_bins:.0f} +/- {np.sqrt(n_data_bins):.1f}")
    print(f"  N_data (entries in mass range):  {n_data_entries}")
    print()
    print("--- signal / background (from fit, in window) ---")
    print(f"  S (Gaussian integral):           {s_fit:.1f} +/- {s_fit_err:.1f}")
    print(f"  B (linear background integral):  {b_fit:.1f}")
    print(f"  S + B (fit):                     {s_fit + b_fit:.1f}  (cf. N_data {n_data_bins:.0f})")
    print(
        f"  yield (full Gaussian, area par): {yield_full:.1f} +/- {yield_full_err:.1f}"
    )
    print(
        f"  peak (predicted): {peak_counts_pred:.1f} counts/bin "
        f"  [density {peak_density:.0f} counts/(GeV/c^2)]"
    )
    print()
    print("--- S/N ---")
    print(f"  S/B:              {sn_ratio:.3f}")
    print(f"  S/sqrt(B):        {significance_b:.2f}")
    print(f"  S/sqrt(S+B):      {significance_tot:.2f}  (S+B = N_data)")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    img_path = os.path.join(script_dir, "../../results/img/lambda_mass_fit.png")
    os.makedirs(os.path.dirname(img_path), exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.stairs(counts, edges, color="k", linewidth=1.0, label="data")
    x_draw = np.linspace(FIT_XMIN, FIT_XMAX, 400)
    y_total = gauss_plus_linear_bin(
        x_draw, bin_width, p["area"], p["center"], p["sigma"], p["slope"], p["intercept"]
    )
    y_gauss = gauss_density(x_draw, p["area"], p["center"], p["sigma"]) * bin_width
    y_linear = linear_density(x_draw, p["center"], p["slope"], p["intercept"]) * bin_width
    ax.plot(x_draw, y_total, "r-", lw=2, label="gauss + linear")
    ax.plot(x_draw, y_gauss, "b--", lw=1.5, label="Gaussian")
    ax.plot(x_draw, y_linear, "g--", lw=1.5, label="linear")
    ax.axvspan(FIT_XMIN, FIT_XMAX, color="r", alpha=0.08)
    ax.set_xlabel(r"$M(p\pi^-)$ [GeV/$c^2$]")
    ax.set_ylabel("Counts / bin")
    ax.set_xlim(1.07, 1.23)
    ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig(img_path, dpi=200, bbox_inches="tight")
    print(f"\nsaved: {img_path}")
    plt.show()


if __name__ == "__main__":
    main()
