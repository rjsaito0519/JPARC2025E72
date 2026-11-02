"""KVC one-photon gain (OPG) per (segment, channel) from the LED-run fits.

Input : OUTPUT_DIR/root/kvc_opg.root (or the path given as the first argument), made by
        analysis/src/kvc_opg.cpp. One entry = one fit of one LED run for (seg, ch), where
          seg     : KVC segment 0-7 (one LED run measures seg s and s+4 together)
          ch      : readout channel 0-3 = a, b, c, d
          board_x : MPPC setting number 1-16 of that channel; settings 4*(seg%4)+1 .. 4*(seg%4)+4
                    are the four MPPCs of the segment.
Output: HDPRM gain lines "6 0 seg 0 ch p0 p1" (p1 = OPG; p0 is a placeholder, the pedestal comes
        from update_hdprm.py) and a plot of the OPG of every MPPC (skip with --no-plot).
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import uproot

project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))
from lib import config
import opg_tool

HV = 58
N_SEG = 8
N_CH = 4
N_MPPC = 4  # MPPCs per (seg, ch)
CH_NAME = {0: "a", 1: "b", 2: "c", 3: "d"}
RESULT_GAIN = 3  # index of the OPG in result_val / result_err


def mppc_index(seg, setting):
    """MPPC 0-3 inside the segment from the setting number 1-16; error if it does not belong to seg."""
    mppc = setting - 1 - N_MPPC * (seg % 4)
    if not 0 <= mppc < N_MPPC:
        raise ValueError(f"setting {setting} does not belong to seg {seg}")
    return mppc


def read_fits(path):
    """{(seg, ch, mppc): [(opg, err), ...]} for the fits taken at HV."""
    tree = uproot.open(path)["tree"].arrays(library="np")
    fits = {}
    for i in range(len(tree["run_num"])):
        if int(tree["hv"][i]) != HV:
            continue
        seg, ch = int(tree["seg"][i]), int(tree["ch"][i])
        setting = int(tree[f"board_{CH_NAME[ch]}"][i])
        key = (seg, ch, mppc_index(seg, setting))
        fits.setdefault(key, []).append((tree["result_val"][i][RESULT_GAIN], tree["result_err"][i][RESULT_GAIN]))
    return fits


def mean_per_mppc(fits):
    """OPG array [seg][ch][mppc][val, err]; weighted mean over repeated fits, NaN if not measured."""
    opg = np.full((N_SEG, N_CH, N_MPPC, 2), np.nan)
    for (seg, ch, mppc), values in fits.items():
        v = np.array(values)
        if len(v) > 1 and opg_tool.check_condition(v[:, 0]):
            print(f"outlier: seg {seg} ch {CH_NAME[ch]} mppc {mppc}: {np.round(v[:, 0], 3)}")
        opg[seg, ch, mppc] = opg_tool.weighted_mean_and_error(v[:, 0], 1.0 / v[:, 1] ** 2)
    return opg


def mean_per_channel(opg):
    """OPG array [seg][ch][val, err]; weighted mean over the four MPPCs of (seg, ch)."""
    gain = np.full((N_SEG, N_CH, 2), np.nan)
    for seg in range(N_SEG):
        for ch in range(N_CH):
            d = opg[seg, ch]
            gain[seg, ch] = opg_tool.weighted_mean_and_error(d[:, 0], 1.0 / d[:, 1] ** 2)
    return gain


def print_missing(opg):
    for seg, ch, mppc in zip(*np.where(np.isnan(opg[:, :, :, 0]))):
        print(f"not measured: seg {seg} ch {CH_NAME[ch]} mppc {mppc}")


def print_hdprm_lines(gain):
    # CId=6 (KVC), PlId=0, SegId=seg, AorT=0 (ADC), UorD=ch, p0=placeholder, p1=OPG
    for seg in range(N_SEG):
        for ch in range(N_CH):
            print(f"6\t0\t{seg}\t0\t{ch}\t{0.0:.6e}\t{gain[seg, ch, 0]:.6e}")


def plot(opg, gain, img_dir):
    fig = plt.figure(figsize=(15, 6))
    ax = fig.add_subplot(111)
    for ch in range(N_CH):
        for seg in range(N_SEG):
            x = N_MPPC * seg + np.arange(1, N_MPPC + 1)
            d = opg[seg, ch]
            ax.errorbar(x, d[:, 0], yerr=d[:, 1], fmt="s", capsize=0, markeredgecolor="k", ms=8,
                        ecolor="black", color=f"C{ch}", markeredgewidth=0.2, zorder=3)
            mean, err = gain[seg, ch]
            ax.hlines(mean, x[0] - 0.5, x[-1] + 0.5, color=f"C{ch}", linestyle="--")
            ax.fill_between([x[0] - 0.5, x[-1] + 0.5], mean - err, mean + err, color=f"C{ch}", alpha=0.2)
        ax.plot([], [], color=f"C{ch}", linestyle="--", label=f"ch {CH_NAME[ch]}")
    ax.set_xticks(N_MPPC * np.arange(N_SEG) + 2.5)
    ax.set_xticklabels([f"seg {s}" for s in range(N_SEG)])
    ax.set_xlabel("segment (4 MPPCs each)")
    ax.set_ylabel("One Photon Gain [arb. unit]")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0))
    plt.subplots_adjust(left=0.1, right=0.85, top=0.98, bottom=0.15)
    img_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(img_dir / f"kvc_opg_{HV}.png", format="png", bbox_inches="tight", dpi=600, transparent=True)
    plt.savefig(img_dir / f"kvc_opg_{HV}.jpg", format="jpg", bbox_inches="tight", dpi=600, transparent=True)


def main():
    for k, v in config.PLOT_SETTINGS.items():
        plt.rcParams[k] = v
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path = Path(args[0]) if args else config.OUTPUT_DIR / "root/kvc_opg.root"
    opg = mean_per_mppc(read_fits(path))
    print_missing(opg)
    gain = mean_per_channel(opg)
    print_hdprm_lines(gain)
    if "--no-plot" not in sys.argv:
        plot(opg, gain, config.OUTPUT_DIR / "img")
        plt.show()


if __name__ == "__main__":
    main()
