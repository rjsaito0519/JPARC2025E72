#!/usr/bin/env python3
"""
run_tpc: TPC の較正のフロントエンド

使い方:
  ./run_tpc.py <root_file> phase  [--mode hit|trk] [--rebin N] [--ofs-update] [--debug]
  ./run_tpc.py <root_file> offset [--mode hit|trk] [--threshold N] [--ave] [--debug]
  ./run_tpc.py <root_file> gain   [--threshold N] [--target-mpv M] [--rebin N] [--de-source all|pion]
                                  [--replace] [--ave] [--debug]
  ./run_tpc.py <helix_root> gain --tree [--threshold N] [--ave] [--debug]
  ./run_tpc.py <root_file> drift  [--mode hit|trk] [--threshold N] [--debug]
  ./run_tpc.py <helix_root> pos   [--model tan|xz] [--threshold N] [--ave] [--parent MAP] [--overwrite] [--debug]

  サブコマンドごとに使えるオプションだけを受け付ける（ほかはエラー）。run 番号は入力 ROOT の
  tpc/run_number の最小値（複数 run の hadd ならいちばん小さい run。読めなければファイル名の runNNNNN）。
  --debug は解析と QA だけ行い param を書かない。

  入力 ROOT の命名（dst_create_runlist / create_runlist と揃える）:
  - phase : hit は run{NNNNN}_TPCHitBcOut.root、trk は TPCTrk_ResY_vs_ClockTime_* を含む DST
  - offset / drift : hit は run{NNNNN}_TPCHitBcOut.root、trk は TPCTrk_ResY_* を含む DST
            （例 run{NNNNN}_TPCHelix.root）。複数 run の hadd も可（いちばん小さい run の TPCPRM を更新）
  - gain : TPCCl_dE_* を含む DST。--tree は π ビームの run の Helix の tpc tree（run{NNNNN}_TPCHelix.root か hadd）
  - pos : run{NNNNN}_TPCHelix.root（inclusive helix DST）または hadd

  - phase : tpc_phase_from_tpcbcout + tpc_phase_plot → TPCPRM の kCobo（aty=3）行 (A, t, w) を
            TpcPhase_CoboFallback から更新。analyzer が使うのはこの行だけで、
            param/TPCPHASE/TpcPhase_NNNNN.root は fit の記録と QA 用
            （--debug では results/root/runNNNNN/runNNNNN_TpcPhase_debug.root に出す）
  - offset: tpc_time_offset_calib（パッドごとの時間オフセット、aty=2 の p0）
  - gain  : tpc_gain_calib（TPCCl_dE_* の Landau fit の MPV → aty=0 の gain）
  - drift : tpc_drift_calib（layer ごとの ResY vs Y 傾き → δv を aty=2 の p1 に加算）
  - pos   : tpc_pos_calib（pad 粒度 TPCPOS。仕様は docs/tpc/TPCPOS_ITER_PREP.ja.md）

  --ave（offset・gain）: fit で更新できなかったパッドを、同じ AsAd の両隣（row±1）の平均、
  両隣が揃わなければ同じ layer・AsAd の更新済みパッドの中央値で埋める（中心の枠は対象外）。
  ただし offset --mode hit では、p0 ではなく fit したパッドの Δp0 を同じ規則で求め、
  今の p0 に足す（BcOut で決まるのはビームの通るパッドだけなので、ほかのパッドは
  trk で決めたパッドごとの差を残したまま全体だけずらす）。
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


def update_tpcprm_phase_stamp(run_num: int, phase_path: Path) -> None:
    """
    TPCPRM の run 別ファイルを更新する（旧 TPCPHASE コメントブロックは除去）。
    ファイルが無ければ TPCParam_example をコピーしてから更新。
    TpcPhase に TpcPhase_CoboFallback があるとき、1 段 fit の kCobo 行を更新。
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
    if phase_path.is_file():
        body = patch_tpcprm_kcobo_from_phase(body, phase_path).rstrip("\n")
    new_text = (body + "\n") if body else ""
    out_prm.write_text(new_text, encoding="utf-8")
    print(colored(f"[INFO] TPCPRM updated: {out_prm}", "green"))


def get_run_number_from_root(root_path: Path) -> int:
    """Smallest tpc/run_number in the file (a hadd of several runs is labelled by its first run).
    Falls back to runNNNNN in the file name."""
    path = Path(root_path).resolve()
    if not path.exists():
        return -1

    try:
        with uproot.open(path) as f:
            if "tpc" in f and "run_number" in f["tpc"]:
                arr = f["tpc"]["run_number"].array(library="np")
                if len(arr) > 0:
                    return int(arr.min())
    except Exception:
        pass

    # Fallback: filename runNNNNN (e.g. run03773_TPCHitBcOut.root / run03773_TPCTracking.root)
    match = re.search(r"run(\d+)", path.name, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return -1


def require_file(path: Path) -> Path:
    path = path.resolve()
    if not path.exists():
        print(colored(f"[Error] File not found: {path}", "red"))
        sys.exit(1)
    return path


def require_binary(bin_dir: Path, name: str) -> Path:
    exe = bin_dir / name
    if not exe.exists():
        print(colored(f"[Error] Binary not found: {exe}. Build myanalysis (./build.sh).", "red"))
        sys.exit(1)
    return exe


def require_run_number(root_path: Path) -> int:
    run_num = get_run_number_from_root(root_path)
    if run_num < 0:
        print(colored("[Error] Could not get run number from ROOT (tpc/run_number) or filename (runNNNNN).", "red"))
        sys.exit(1)
    return run_num


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="TPC calibration frontend. Usage: run_tpc.py <root> {phase,offset,gain,drift,pos} [options]"
    )
    parser.add_argument(
        "root_file",
        type=Path,
        help="Input ROOT: phase=runNNNNN_TPCHitBcOut.root; offset/drift=TPCHitBcOut (hit) or a DST with "
        "TPCTrk_ResY_* such as runNNNNN_TPCHelix.root (trk); gain=DST with TPCCl_dE_*; "
        "gain --tree / pos=runNNNNN_TPCHelix.root.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True, metavar="{phase,offset,gain,drift,pos}")

    def add_debug(p: argparse.ArgumentParser) -> None:
        p.add_argument("--debug", action="store_true", help="Run the analysis and QA but do not write params.")

    def add_mode(p: argparse.ArgumentParser, text: str) -> None:
        p.add_argument("--mode", choices=["hit", "trk"], default="hit", help=text)

    p = sub.add_parser("phase", help="TPC clock phase (tpc_phase_from_tpcbcout + tpc_phase_plot)")
    add_mode(p, "Histogram preference (falls back to the other if missing).")
    p.add_argument("--rebin", type=int, default=None, metavar="N",
                   help="Rebin the 2D clock axis before the gaus profile.")
    p.add_argument("--ofs-update", action="store_true",
                   help="Apply the baseline c0 as a CoBo-uniform time offset (dp0 = -c0; skips the center-frame X). "
                   "Ignored with --debug.")
    add_debug(p)

    p = sub.add_parser("offset", help="Per-pad time offset (tpc_time_offset_calib)")
    add_mode(p, "hit / trk histogram set.")
    p.add_argument("--threshold", type=int, default=None, metavar="N",
                   help="Minimum entries in the fit histogram of a pad.")
    p.add_argument("--ave", action="store_true",
                   help="Fill pads without a fit; with --mode hit, shift them by the fitted pads' dp0 "
                        "(see the module help).")
    p.add_argument("--common", action="store_true",
                   help="--mode hit only, instead of --ave: shift every pad of an AsAd by the median dp0 of its "
                        "fitted pads (global median for AsAds with too few), keeping the pad-to-pad differences.")
    p.add_argument("--ywin", type=float, default=None, metavar="Y",
                   help="--mode trk only: fit only the hits with |y| < Y mm (the height of the BcOut beam).")
    p.add_argument("--ywin-max-layer", type=int, default=None, metavar="L",
                   help="Apply --ywin only to layers <= L.")
    add_debug(p)

    p = sub.add_parser("gain", help="Per-pad gain (tpc_gain_calib)")
    p.add_argument("--threshold", type=int, default=None, metavar="N",
                   help="Minimum entries in the dE histogram of a pad.")
    p.add_argument("--target-mpv", type=float, default=None, metavar="M",
                   help="Target MPV (default 200). Histogram mode only.")
    p.add_argument("--rebin", type=int, default=None, metavar="N", help="Rebin of the dE histogram. Histogram mode only.")
    p.add_argument("--de-source", choices=["all", "pion"], default=None,
                   help="dE histograms: all (TPCCl_dE_*) or pion (TPCCl_dE_Pion_*, default). Histogram mode only.")
    p.add_argument("--tree", action="store_true",
                   help="Pion sample from the Helix tpc tree of pion-beam runs (see common/include/tpc_gain_tree.h).")
    p.add_argument("--replace", action="store_true",
                   help="Set gain = target/MPV instead of multiplying the previous gain. Histogram mode only.")
    p.add_argument("--ave", action="store_true", help="Fill pads without a fit (see the module help).")
    add_debug(p)

    p = sub.add_parser("drift", help="Drift velocity per layer (tpc_drift_calib)")
    add_mode(p, "Histogram preference (falls back to the other if missing).")
    p.add_argument("--threshold", type=int, default=None, metavar="N",
                   help="Minimum entries of a layer 2D histogram (default 100000).")
    add_debug(p)

    p = sub.add_parser("pos", help="TPCPOS pad-level fit (tpc_pos_calib)")
    p.add_argument("--model", default=None, metavar="M", help="Update model: tan (default) or xz.")
    p.add_argument("--threshold", type=int, default=None, metavar="N", help="Minimum hits per pad (default 100).")
    p.add_argument("--ave", action="store_true",
                   help="Fill unfitted pads (layer+ASAD mean; center-frame forced to (0,0)).")
    p.add_argument("--parent", type=Path, default=None, metavar="MAP",
                   help="TPCPOS map used to reconstruct the helix input; the update is added onto it.")
    p.add_argument("--overwrite", action="store_true",
                   help="Replace touched pads by -<r> instead of adding to the parent map.")
    add_debug(p)
    return parser


def opt(flag: str, value) -> list[str]:
    return [] if value is None else [flag, str(value)]


def run_phase(args, bin_dir: Path) -> None:
    bcout_path = require_file(args.root_file)
    run_num = require_run_number(bcout_path)
    print(colored(f"[INFO] Run number: {run_num}", "green"))
    print(colored(f"[INFO] Input file: {bcout_path}", "green"))
    exe_fit = require_binary(bin_dir, "tpc_phase_from_tpcbcout")
    exe_plot = require_binary(bin_dir, "tpc_phase_plot")

    # --debug leaves param/TPCPHASE untouched and writes the fit result next to the other outputs.
    if args.debug:
        phase_dir = Path(config.OUTPUT_DIR) / "root" / f"run{run_num:05d}"
        phase_path = phase_dir / f"run{run_num:05d}_TpcPhase_debug.root"
    else:
        phase_dir = config.PARAM_DIR / "TPCPHASE"
        phase_path = phase_dir / f"TpcPhase_{run_num:05d}.root"
    phase_dir.mkdir(parents=True, exist_ok=True)

    cmd_fit = [str(exe_fit), str(bcout_path), str(phase_path), "--run", str(run_num), "--mode", args.mode,
               *opt("--rebin", args.rebin)]
    if args.ofs_update and not args.debug:
        cmd_fit.append("--ofs-update")
    print(colored(">>> Step 1: TPC phase fit (tpc_phase_from_tpcbcout)", "cyan"))
    run_command(" ".join(cmd_fit))

    pdf_path = Path(config.OUTPUT_DIR) / "img" / f"run{run_num:05d}" / f"{bcout_path.stem}_TPC_Phase.pdf"
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    cmd_plot = [str(exe_plot), str(bcout_path), "--phase", str(phase_path), "--mode", args.mode, str(pdf_path)]
    print(colored(">>> Step 2: TPC phase plot (tpc_phase_plot)", "cyan"))
    run_command(" ".join(cmd_plot))

    if not args.debug:
        update_tpcprm_phase_stamp(run_num, phase_path)

    print(colored(f"\n[DONE] TPC phase calibration complete for run {run_num}.", "green", attrs=["bold"]))
    print(colored(f"  TpcPhase: {phase_path}", "green"))
    print(colored(f"  PDF: {pdf_path}", "green"))


def run_offset(args, bin_dir: Path) -> None:
    input_root = require_file(args.root_file)
    run_num = require_run_number(input_root)
    exe = require_binary(bin_dir, "tpc_time_offset_calib")
    cmd = [str(exe), str(input_root), str(run_num), "--mode", args.mode, *opt("--threshold", args.threshold)]
    cmd += ["--ave"] if args.ave else []
    cmd += ["--common"] if args.common else []
    cmd += [*opt("--ywin", args.ywin), *opt("--ywin-max-layer", args.ywin_max_layer)]
    cmd += ["--debug"] if args.debug else []
    print(colored(f">>> TPC time offset calib (tpc_time_offset_calib), run={run_num}", "cyan"))
    run_command(" ".join(cmd))


def run_gain(args, bin_dir: Path) -> None:
    input_root = require_file(args.root_file)
    run_num = require_run_number(input_root)
    exe = require_binary(bin_dir, "tpc_gain_calib")
    cmd = [str(exe), str(input_root), str(run_num), *opt("--de-source", args.de_source),
           *opt("--threshold", args.threshold), *opt("--target-mpv", args.target_mpv), *opt("--rebin", args.rebin)]
    for flag in ("tree", "replace", "ave", "debug"):
        if getattr(args, flag):
            cmd.append(f"--{flag}")
    print(colored(f">>> TPC gain calib (tpc_gain_calib), run={run_num}", "cyan"))
    run_command(" ".join(cmd))


def run_drift(args, bin_dir: Path) -> None:
    input_root = require_file(args.root_file)
    run_num = require_run_number(input_root)
    exe = require_binary(bin_dir, "tpc_drift_calib")
    cmd = [str(exe), str(input_root), str(run_num), "--mode", args.mode, *opt("--threshold", args.threshold)]
    cmd += ["--debug"] if args.debug else []
    print(colored(f">>> TPC drift calib (tpc_drift_calib), run={run_num}", "cyan"))
    run_command(" ".join(cmd))


def run_pos(args, bin_dir: Path) -> None:
    # Spec: docs/tpc/TPCPOS_ITER_PREP.ja.md
    input_root = require_file(args.root_file)
    exe = require_binary(bin_dir, "tpc_pos_calib")
    cmd = [str(exe), str(input_root), *opt("--model", args.model), *opt("--threshold", args.threshold)]
    cmd += ["--ave"] if args.ave else []
    cmd += ["--debug"] if args.debug else []
    if args.parent is not None:
        parent_map = args.parent.resolve()
        if not parent_map.exists():
            print(colored(f"[Error] parent map not found: {parent_map}", "red"))
            sys.exit(1)
        cmd += ["--parent", str(parent_map)]
    cmd += ["--overwrite"] if args.overwrite else []
    print(colored(">>> TPC pos calib (tpc_pos_calib)", "cyan"))
    run_command(" ".join(cmd))


def main():
    args = build_parser().parse_args()
    bin_dir = project_root / "bin"
    {"phase": run_phase, "offset": run_offset, "gain": run_gain, "drift": run_drift, "pos": run_pos}[args.cmd](
        args, bin_dir
    )


if __name__ == "__main__":
    main()

