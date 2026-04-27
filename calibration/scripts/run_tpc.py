#!/usr/bin/env python3
"""
run_tpc: TPC 周りの便利フロントエンド

使い方:
  ./run_tpc.py <root_file> phase  [--fit-step N] [--vdrift V] [--smooth N] [--mode hit|trk]
  ./run_tpc.py <root_file> phase --plot-only [--vdrift V]
  ./run_tpc.py <root_file> offset [--run N] [--mode hit|trk] [...]
  ./run_tpc.py <root_file> gain [--run N] [--target-mpv 200] [...]

  - phase : tpc_phase_from_tpcbcout + tpc_phase_plot（root は tpcbcout）
            → --fit-step 1 なら kCobo(3) 行を TTree から更新（TPCPRM へのコメント追記はしない）
  - offset: tpc_time_offset_calib（root は tpc_runXXXXX 等）
            → QA PDF に p0 / Δp0 の分布（1D + TPC 上の 2D マップ）を含む
  - gain  : tpc_gain_calib（root は tpc_runXXXXX 等）
            → TPCCl_dE_* の Landau fit から MPV を求めて ATY=0 の gain を更新
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import uproot
from termcolor import colored

# Add shared library path to find config
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config


def run_command(cmd: str, exit_on_error: bool = True) -> int:
    print(colored(f"Running: {cmd}", "cyan"))
    ret = subprocess.call(cmd, shell=True)
    if ret != 0:
        print(colored(f"[Error] Command failed with exit code {ret}", "red"))
        if exit_on_error:
            sys.exit(ret)
    return ret


# 互換のため旧スタンプ文字列は保持（更新時に除去する）
TPCPHASE_STAMP_BEGIN = "# --- begin TPCPHASE pointer (run_tpc.py; not parsed by TPCParamMan) ---"
TPCPHASE_STAMP_END = "# --- end TPCPHASE pointer ---"


KCOBO_LINE_PREFIX = re.compile(r"^(\d+)[\t ]+0[\t ]+3[\t ]+")


def patch_tpcprm_kcobo_from_phase(text: str, phase_path: Path) -> str:
    """
    TpcPhase_*.root の TTree TpcPhase_CoboFallback から、--fit-step 1 の CoBo だけ
    TPCPRM の「cobo 0 3 p0 p1 p2」行を TPCParamMan::PhaseShift フォールバックと整合するよう更新。
    nstep_fit==0（flat 初期）→ 0,0,1。多段や失敗（nstep_fit!=1）の CoBo は行を触らない。
    """
    updates: dict[int, tuple[float, float, float]] = {}
    try:
        with uproot.open(phase_path) as pf:
            tree_key = None
            for k in pf.keys():
                if str(k).split(";")[0] == "TpcPhase_CoboFallback":
                    tree_key = k
                    break
            if tree_key is None:
                return text
            t = pf[tree_key]
            d = t.arrays(["cobo", "nstep_fit", "p0", "p1", "p2"], library="np")
    except Exception as exc:
        print(
            colored(
                f"[WARN] kCobo lines not patched (read TpcPhase_CoboFallback failed): {exc}",
                "yellow",
            )
        )
        return text

    n = len(d["cobo"])
    for i in range(n):
        c = int(d["cobo"][i])
        ns = int(d["nstep_fit"][i])
        if ns == 1:
            updates[c] = (float(d["p0"][i]), float(d["p1"][i]), float(d["p2"][i]))
        elif ns == 0:
            updates[c] = (0.0, 0.0, 1.0)
    if not updates:
        return text

    lines = text.splitlines(keepends=True)
    out: list[str] = []
    n_patched = 0
    for line in lines:
        m = KCOBO_LINE_PREFIX.match(line)
        if m:
            c = int(m.group(1))
            if c in updates:
                p0, p1, p2 = updates[c]
                nl = "\n"
                if line.endswith("\r\n"):
                    nl = "\r\n"
                elif line.endswith("\r"):
                    nl = "\r"
                out.append(f"{c}\t0\t3\t{p0:.8f}\t{p1:.8f}\t{p2:.8f}{nl}")
                n_patched += 1
                continue
        out.append(line)
    print(colored(f"[INFO] TPCPRM kCobo (aty=3) lines patched: {n_patched} (from phase tree)", "green"))
    return "".join(out)


def phase_path_for_analyzer_conf(phase_path: Path, param_dir: Path) -> str:
    """analyzer conf の TPCPHASE: にそのまま書ける形式（param/TPCPHASE/...）。"""
    phase_res = phase_path.resolve()
    param_res = param_dir.resolve()
    try:
        rel = phase_res.relative_to(param_res)
        return f"param/{rel.as_posix()}"
    except ValueError:
        return str(phase_res)


def update_tpcprm_phase_stamp(
    run_num: int, phase_path: Path, *, patch_kcobo: bool = True
) -> None:
    """
    TPCPRM の run 別ファイルを更新する（旧 TPCPHASE コメントブロックは除去）。
    ファイルが無ければ TPCParam_0_yoffset_adjusted をコピーしてから更新。
    patch_kcobo=True かつ TpcPhase に TpcPhase_CoboFallback があるとき、--fit-step 1 の kCobo 行を更新。
    """
    param_dir = config.PARAM_DIR
    tpcprm_e72 = param_dir / "TPCPRM" / "e72"
    tpcprm_e72.mkdir(parents=True, exist_ok=True)
    out_prm = tpcprm_e72 / f"TPCParam_e72_run{run_num:05d}"
    base_prm = param_dir / "TPCPRM" / "TPCParam_example"


    if not out_prm.exists():
        if not base_prm.is_file():
            print(
                colored(
                    f"[WARN] Skip TPCPRM update: no {out_prm.name} and missing base {base_prm}",
                    "yellow",
                )
            )
            return
        shutil.copy2(base_prm, out_prm)
        print(colored(f"[INFO] Created {out_prm.name} from {base_prm.name}", "green"))

    raw = out_prm.read_text(encoding="utf-8", errors="replace")
    lines = raw.splitlines(keepends=True)
    out_lines: list[str] = []
    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if stripped == TPCPHASE_STAMP_BEGIN.strip():
            i += 1
            while i < len(lines) and lines[i].strip() != TPCPHASE_STAMP_END.strip():
                i += 1
            if i < len(lines):
                i += 1
            continue
        out_lines.append(lines[i])
        i += 1

    body = "".join(out_lines).rstrip("\n")
    if patch_kcobo and phase_path.is_file():
        body = patch_tpcprm_kcobo_from_phase(body, phase_path).rstrip("\n")
    new_text = (body + "\n") if body else ""
    out_prm.write_text(new_text, encoding="utf-8")
    print(colored(f"[INFO] TPCPRM updated: {out_prm}", "green"))


def get_run_number_from_root(root_path: Path) -> int:
    """ROOT ファイルの tpc ツリー run_number から取得。失敗時はファイル名の runNNNNN にフォールバック。"""
    path = Path(root_path).resolve()
    if not path.exists():
        return -1

    try:
        with uproot.open(path) as f:
            if "tpc" not in f:
                return -1
            tree = f["tpc"]
            if "run_number" not in tree:
                return -1
            arr = tree["run_number"].array(library="np")
            if arr is None or len(arr) == 0:
                return -1
            return int(arr[0])
    except Exception:
        pass

    # Fallback: filename runNNNNN
    match = re.search(r"run(\d+)", path.name, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return -1


def main():
    parser = argparse.ArgumentParser(
        description="TPC helper frontend. Usage: run_tpc.py <root> phase|offset|gain [options]"
    )
    parser.add_argument(
        "root_file",
        type=Path,
        help="Input ROOT file (tpcbcout for phase, tpc_runXXXXX for offset).",
    )
    parser.add_argument(
        "cmd",
        choices=["phase", "offset", "gain"],
        help="Subcommand: phase (tpc_phase_from_tpcbcout + plot), offset (tpc_time_offset_calib), or gain (tpc_gain_calib).",
    )
    # phase 用オプション
    parser.add_argument(
        "--fit-step",
        type=int,
        default=1,
        metavar="N",
        help="[phase] Step function order for fit (default: 1). Use 0 for zero-initial TpcPhase only.",
    )
    parser.add_argument(
        "--vdrift",
        type=float,
        default=None,
        metavar="V",
        help="[phase] Drift velocity [mm/ns]. Passed to both programs if set.",
    )
    parser.add_argument(
        "--smooth",
        type=int,
        default=None,
        metavar="N",
        help="[phase] Moving average half-window for tpc_phase_from_tpcbcout.",
    )
    parser.add_argument(
        "--step-width",
        type=float,
        default=None,
        metavar="W",
        help="[phase] Fixed step width [ns] for fit function (smaller => steeper edge).",
    )
    parser.add_argument(
        "--free",
        action="store_true",
        help="[phase] Free the step width in fit (default: fixed).",
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="[phase] Skip Step1 (tpc_phase_from_tpcbcout) and run plot-only using existing TpcPhase_*.root.",
    )
    parser.add_argument(
        "--phase-root",
        type=Path,
        default=None,
        metavar="FILE",
        help="[phase] Override TpcPhase root file path (used with --plot-only).",
    )
    prm_group = parser.add_mutually_exclusive_group()
    prm_group.add_argument(
        "--update-tpcprm",
        dest="update_tpcprm",
        action="store_true",
        default=None,
        help="[phase] Stamp TPCPHASE path into TPCPRM/e72/TPCParam_e72_runNNNNN (default: on unless --plot-only).",
    )
    prm_group.add_argument(
        "--no-update-tpcprm",
        dest="update_tpcprm",
        action="store_false",
        default=None,
        help="[phase] Do not modify TPCPRM.",
    )
    parser.add_argument(
        "--no-patch-kcobo",
        action="store_true",
        help="[phase] Do not rewrite kCobo (aty=3) lines from TpcPhase_CoboFallback (TPCPRM stamp only).",
    )
    # offset 用オプション
    parser.add_argument(
        "--run",
        type=int,
        default=None,
        help="[offset] Run number (if omitted, try to detect from ROOT file).",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["hit", "trk"],
        default="hit",
        help="hit / trk: [offset] TPCHit vs TPCTrk ResY hist; [phase] ClockTime vs ResY 2D の名前優先順。",
    )
    parser.add_argument(
        "--threshold",
        type=int,
        default=None,
        help="[offset] Minimum entries in ±stat-range window (passed to C++).",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=None,
        help="[offset] Upper limit for fitted sigma (passed to C++ --sigma).",
    )
    parser.add_argument(
        "--nsigma",
        type=float,
        default=None,
        help="[offset] Fit range = peak ± nsigma * RMS (passed to C++ --nsigma).",
    )
    parser.add_argument(
        "--stat-range",
        type=float,
        default=None,
        help="[offset] Statistics range around 0 for entry counting (passed to C++).",
    )
    parser.add_argument(
        "--min-peak-ratio",
        type=float,
        default=None,
        metavar="R",
        help="[offset] Skip fit if peak < R * mean in window (flat hist rejection, passed to C++ --min-peak-ratio).",
    )
    parser.add_argument(
        "--min-sigma",
        type=float,
        default=None,
        metavar="MM",
        help="[offset] Reject fit if Gaussian sigma < this (mm); also used for SetParLimits floor (C++ --min-sigma).",
    )
    parser.add_argument(
        "--min-sigma-rel-rms",
        type=float,
        default=None,
        metavar="F",
        help="[offset] Also require sigma >= F * hist RMS (needle rejection; C++ --min-sigma-rel-rms).",
    )
    parser.add_argument(
        "--sigma-init-rel",
        type=float,
        default=None,
        metavar="R",
        help="[offset] Initial Gaussian sigma = R * max(h_rms, hist RMS) (C++ --sigma-init-rel, default ~0.72).",
    )
    parser.add_argument(
        "--max-chi2-ndf",
        type=float,
        default=None,
        metavar="X",
        help="[offset] Reject fit if chi2/ndf > X (0 disables; C++ --max-chi2-ndf).",
    )
    parser.add_argument(
        "--min-fit-ndf",
        type=int,
        default=None,
        help="[offset] Require fit NDF >= N (C++ --min-fit-ndf).",
    )
    parser.add_argument(
        "--min-fit-prob",
        type=float,
        default=None,
        metavar="P",
        help="[offset] Require ROOT TMath::Prob(chi2,ndf) >= P (0 disables; C++ --min-fit-prob).",
    )
    parser.add_argument(
        "--local-peak-mm",
        type=float,
        default=None,
        metavar="W",
        help="[offset] Half-width [mm] for local RMS around tallest peak (0=use full hist RMS; C++ --local-peak-mm).",
    )
    parser.add_argument(
        "--local-peak-sep",
        type=int,
        default=None,
        metavar="B",
        help="[offset] Peak locality in bins (±B); taller peak wins, ties→|x| smaller (C++ --local-peak-sep).",
    )
    parser.add_argument(
        "--ave",
        action="store_true",
        help="[offset] Fill non-fitted pads with mean p0 of fitted pads in the same ASAD (exclude center-frame pads).",
    )
    # gain 用オプション
    parser.add_argument(
        "--target-mpv",
        type=float,
        default=None,
        metavar="M",
        help="[gain] Target MPV for gain scaling (C++ --target-mpv, default 200.0).",
    )
    parser.add_argument(
        "--fit-nsigma",
        type=float,
        default=None,
        metavar="N",
        help="[gain] Fit range = peak ± nsigma * local RMS (C++ --fit-nsigma).",
    )
    parser.add_argument(
        "--local-peak-half",
        type=float,
        default=None,
        metavar="W",
        help="[gain] Half width [ADC] for local RMS around peak (C++ --local-peak-half).",
    )
    parser.add_argument(
        "--min-width",
        type=float,
        default=None,
        metavar="W",
        help="[gain] Minimum Landau width accepted (C++ --min-width).",
    )
    parser.add_argument(
        "--rebin",
        type=int,
        default=None,
        metavar="N",
        help="[gain] Rebin factor for TPCCl_dE histograms before fit (C++ --rebin).",
    )
    parser.add_argument(
        "--mpv-min",
        type=float,
        default=None,
        metavar="X",
        help="[gain] Minimum allowed MPV for fit acceptance/seed (C++ --mpv-min, default 100).",
    )
    parser.add_argument(
        "--mpv-max",
        type=float,
        default=None,
        metavar="X",
        help="[gain] Maximum allowed MPV for fit acceptance/seed (C++ --mpv-max; <=0 disables).",
    )
    de_source_group = parser.add_mutually_exclusive_group()
    de_source_group.add_argument(
        "--de-source",
        type=str,
        choices=["all", "pion"],
        default=None,
        help="[gain] dE histogram source: all (TPCCl_dE_*) or pion (TPCCl_dE_Pion_*).",
    )
    de_source_group.add_argument(
        "--pion",
        action="store_true",
        help="[gain] Shortcut of --de-source pion.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Dry-run mode: run fit/plot but do not update parameter files.",
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="[gain] Replace gain by target_mpv/MPV instead of multiplying old gain.",
    )

    parser.set_defaults(update_tpcprm=None)
    args = parser.parse_args()

    bin_dir = project_root / "bin"

    if args.cmd == "phase":
        bcout_path = args.root_file.resolve()
        if not bcout_path.exists():
            print(colored(f"[Error] File not found: {bcout_path}", "red"))
            sys.exit(1)

        run_num = get_run_number_from_root(bcout_path)
        if run_num < 0:
            print(
                colored(
                    "[Error] Could not get run number from ROOT (tpc/run_number) or filename (runNNNNN).",
                    "red",
                )
            )
            sys.exit(1)
        print(colored(f"[INFO] Run number: {run_num}", "green"))
        print(colored(f"[INFO] Input file: {bcout_path}", "green"))

        phase_dir = config.PARAM_DIR / "TPCPHASE"
        phase_dir.mkdir(parents=True, exist_ok=True)
        phase_path = (
            args.phase_root.resolve()
            if args.phase_root is not None
            else (phase_dir / f"TpcPhase_{run_num:05d}.root")
        )

        exe_fit = bin_dir / "tpc_phase_from_tpcbcout"
        exe_plot = bin_dir / "tpc_phase_plot"
        if not exe_fit.exists():
            print(colored(f"[Error] Binary not found: {exe_fit}. Please build the project.", "red"))
            sys.exit(1)
        if not exe_plot.exists():
            print(colored(f"[Error] Binary not found: {exe_plot}. Please build the project.", "red"))
            sys.exit(1)

        if args.plot_only and not phase_path.exists():
            print(colored(f"[Error] --plot-only requested but phase file not found: {phase_path}", "red"))
            print(colored("Run phase once (without --plot-only) to generate it, or pass --phase-root.", "red"))
            sys.exit(1)

        vdrift_opt = f" --vdrift {args.vdrift}" if args.vdrift is not None else ""
        smooth_opt = f" --smooth {args.smooth}" if args.smooth is not None else ""
        mode_opt = f" --mode {args.mode}"

        if not args.plot_only:
            # Step 1: tpc_phase_from_tpcbcout
            if args.fit_step == 0:
                cmd_fit = f"{exe_fit} {phase_path}{vdrift_opt} --fit-step 0{mode_opt}"
            else:
                cmd_fit = (
                    f"{exe_fit} {bcout_path} {phase_path} --fit-step {args.fit_step}"
                    f"{vdrift_opt}{smooth_opt}{mode_opt}"
                )
            if args.step_width is not None:
                cmd_fit += f" --step-width {args.step_width}"
            if args.free:
                cmd_fit += " --free"
            print(colored(">>> Step 1: TPC phase fit (tpc_phase_from_tpcbcout)", "cyan"))
            run_command(cmd_fit)

        # Step 2: tpc_phase_plot
        cmd_plot = f"{exe_plot} {bcout_path} --phase {phase_path}{vdrift_opt}{mode_opt}"
        print(colored(">>> Step 2: TPC phase plot (tpc_phase_plot)", "cyan"))
        run_command(cmd_plot)

        do_tpcprm = args.update_tpcprm
        if do_tpcprm is None:
            do_tpcprm = not args.plot_only
        if args.debug:
            do_tpcprm = False
        if do_tpcprm:
            if not phase_path.is_file():
                print(
                    colored(
                        f"[WARN] Skip TPCPRM stamp: phase file not found: {phase_path}",
                        "yellow",
                    )
                )
            else:
                update_tpcprm_phase_stamp(
                    run_num, phase_path, patch_kcobo=not args.no_patch_kcobo
                )

        print(
            colored(
                f"\n[DONE] TPC phase calibration complete for run {run_num}.",
                "green",
                attrs=["bold"],
            )
        )
        print(colored(f"  TpcPhase: {phase_path}", "green"))
        print(
            colored(
                f"  PDF: {config.OUTPUT_DIR}/img/run{run_num:05d}/  (see C++ output)",
                "green",
            )
        )

    elif args.cmd == "offset":
        input_root = args.root_file.resolve()
        if not input_root.exists():
            print(colored(f"[Error] File not found: {input_root}", "red"))
            sys.exit(1)

        run_num = args.run if args.run is not None else get_run_number_from_root(input_root)
        if run_num < 0:
            print(
                colored(
                    "[Error] Could not get run number (use --run or ensure tpc/run_number is present).",
                    "red",
                )
            )
            sys.exit(1)

        exe_offset = bin_dir / "tpc_time_offset_calib"
        if not exe_offset.exists():
            print(colored(f"[Error] Binary not found: {exe_offset}. Please build the project.", "red"))
            sys.exit(1)

        opts = [str(exe_offset), str(input_root), str(run_num)]
        if args.mode:
            opts += ["--mode", args.mode]
        if args.threshold is not None:
            opts += ["--threshold", str(args.threshold)]
        if args.sigma is not None:
            opts += ["--sigma", str(args.sigma)]
        if args.nsigma is not None:
            opts += ["--nsigma", str(args.nsigma)]
        if args.stat_range is not None:
            opts += ["--stat-range", str(args.stat_range)]
        if args.min_peak_ratio is not None:
            opts += ["--min-peak-ratio", str(args.min_peak_ratio)]
        if args.min_sigma is not None:
            opts += ["--min-sigma", str(args.min_sigma)]
        if args.min_sigma_rel_rms is not None:
            opts += ["--min-sigma-rel-rms", str(args.min_sigma_rel_rms)]
        if args.sigma_init_rel is not None:
            opts += ["--sigma-init-rel", str(args.sigma_init_rel)]
        if args.max_chi2_ndf is not None:
            opts += ["--max-chi2-ndf", str(args.max_chi2_ndf)]
        if args.min_fit_ndf is not None:
            opts += ["--min-fit-ndf", str(args.min_fit_ndf)]
        if args.min_fit_prob is not None:
            opts += ["--min-fit-prob", str(args.min_fit_prob)]
        if args.local_peak_mm is not None:
            opts += ["--local-peak-mm", str(args.local_peak_mm)]
        if args.local_peak_sep is not None:
            opts += ["--local-peak-sep", str(args.local_peak_sep)]
        if args.ave:
            opts += ["--ave"]
        if args.debug:
            opts += ["--debug"]

        cmd = " ".join(opts)
        print(colored(f">>> TPC time offset calib (tpc_time_offset_calib), run={run_num}", "cyan"))
        run_command(cmd)

    elif args.cmd == "gain":
        input_root = args.root_file.resolve()
        if not input_root.exists():
            print(colored(f"[Error] File not found: {input_root}", "red"))
            sys.exit(1)

        run_num = args.run if args.run is not None else get_run_number_from_root(input_root)
        if run_num < 0:
            print(
                colored(
                    "[Error] Could not get run number (use --run or ensure tpc/run_number is present).",
                    "red",
                )
            )
            sys.exit(1)

        exe_gain = bin_dir / "tpc_gain_calib"
        if not exe_gain.exists():
            print(colored(f"[Error] Binary not found: {exe_gain}. Please build the project.", "red"))
            sys.exit(1)

        opts = [str(exe_gain), str(input_root), str(run_num)]
        if args.mode:
            opts += ["--mode", args.mode]
        if args.target_mpv is not None:
            opts += ["--target-mpv", str(args.target_mpv)]
        if args.threshold is not None:
            opts += ["--threshold", str(args.threshold)]
        if args.stat_range is not None:
            opts += ["--stat-range", str(args.stat_range)]
        if args.min_peak_ratio is not None:
            opts += ["--min-peak-ratio", str(args.min_peak_ratio)]
        if args.max_chi2_ndf is not None:
            opts += ["--max-chi2-ndf", str(args.max_chi2_ndf)]
        if args.min_fit_ndf is not None:
            opts += ["--min-fit-ndf", str(args.min_fit_ndf)]
        if args.min_fit_prob is not None:
            opts += ["--min-fit-prob", str(args.min_fit_prob)]
        if args.fit_nsigma is not None:
            opts += ["--fit-nsigma", str(args.fit_nsigma)]
        if args.local_peak_half is not None:
            opts += ["--local-peak-half", str(args.local_peak_half)]
        if args.min_width is not None:
            opts += ["--min-width", str(args.min_width)]
        if args.rebin is not None:
            opts += ["--rebin", str(args.rebin)]
        if args.mpv_min is not None:
            opts += ["--mpv-min", str(args.mpv_min)]
        if args.mpv_max is not None:
            opts += ["--mpv-max", str(args.mpv_max)]
        de_source = "all"
        if args.pion:
            de_source = "pion"
        elif args.de_source is not None:
            de_source = args.de_source
        opts += ["--de-source", de_source]
        if args.debug:
            opts += ["--debug"]
        if args.replace:
            opts += ["--replace"]

        cmd = " ".join(opts)
        print(colored(f">>> TPC gain calib (tpc_gain_calib), run={run_num}", "cyan"))
        run_command(cmd)


if __name__ == "__main__":
    main()

