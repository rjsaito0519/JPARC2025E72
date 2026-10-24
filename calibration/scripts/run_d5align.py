#!/usr/bin/env python3
"""BLC1 alignment against BLC2 through the D5 matrix (BLC2 fixed, BLC1a/b translated).

Modes (run numbers as positional arguments, several runs allowed):

  prep    Make the decode inputs for the BLC2-referenced D5 fit:
            DCGEO  param/DCGEO/e72/DCGeomParam_runXXXXX_K_bc1reso  (BLC1a/b Res x FACTOR)
            conf   param/conf/e72/analyzer_runXXXXX_d5_K_bc1reso.conf
            runlist runmanager/runlist/scratch/d5bc1reso_<runs>.yml (+ d5_<runs>.yml)
          and print the myrun commands. Re-run after every DCGEO update.
  review  (convergence table -> /group/had/sks/Users/sryuta/tmp/d5_align/)
          y: d5_bc1reso_review --<beam> on runXXXXX_D5_<S>_bc1reso.root (BLC1a/b dy, BLC2-referenced)
          x: d5_fitfree_review on runXXXXX_D5_K.root (BLC1a/b opposite shift s)
  apply   Back up DCGEO and write the review result into Ofs with apply_dcgeo_offset.py
          (dry-run unless --write). --y / --x select the components (default both).

One iteration:
  run_d5align.py prep <runs>  ->  myrun d5 + d5bc1reso  ->  run_d5align.py review <runs>
  -> run_d5align.py apply <runs> --write  ->  myrun bcin  ->  run_dcprm.py <run> resi --bcin [--kaon]
  -> (prep again) ...  until |dy| and |s| are below the tolerance.

The common x shift of BLC1 is degenerate with delta and is left unchanged.
"""

import argparse
import datetime
import re
import shutil
import subprocess
import sys
from pathlib import Path
from termcolor import colored

project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config

DCGEO_DIR = config.PARAM_DIR / "DCGEO" / config.SUB_DIR
CONF_DIR = config.PARAM_DIR / "conf" / config.SUB_DIR
RUNLIST_DIR = config.ANALYZER_DIR / "runmanager" / "runlist"
BIN_DIR = project_root / "bin"
APPLY_SCRIPT = Path(__file__).resolve().parent / "apply_dcgeo_offset.py"
BACKUP_BASE = Path("/group/had/sks/Users/sryuta/backup")
BEAM = {"kaon": ("K", "--kaon"), "pion": ("Pi", "--pion")}  # --beam -> (param suffix, review flag)
SUFFIX = "K"  # set from --beam in main()
BEAM_FLAG = "--kaon"
TAG = "bc1reso"
SUMMARY_DIR = Path("/group/had/sks/Users/sryuta/tmp/d5_align")  # per-iteration convergence tables (scratch)


def run_tag(run):
    return f"run{run:05d}"


def dcgeo_path(run):
    return DCGEO_DIR / f"DCGeomParam_{run_tag(run)}_{SUFFIX}"


def d5_root(run, tag=None):
    name = f"{run_tag(run)}_D5_{SUFFIX}" + (f"_{tag}" if tag else "") + ".root"
    return config.DECODE_DIR / run_tag(run) / name


def txt_dir(run):
    return SUMMARY_DIR / run_tag(run)


def run_cmd(cmd, capture=False):
    print(colored(f"Running: {cmd}", "cyan"))
    if capture:
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return res.returncode, res.stdout + res.stderr
    return subprocess.call(cmd, shell=True), ""


def runs_label(runs):
    return f"{runs[0]:05d}-{runs[-1]:05d}_{len(runs)}runs" if len(runs) > 1 else f"{runs[0]:05d}"


# -- prep ---------------------------------------------------------------------
def inflate_blc1_res(src, dst, factor):
    out = []
    for line in src.read_text().splitlines(keepends=True):
        f = line.split()
        if len(f) >= 13 and f[0].isdigit() and re.match(r"BLC1[ab]-", f[1]):
            f[9] = f"{float(f[9]) * factor:.6g}"
            end = "\n" if line.endswith("\n") else ""
            line = "\t".join(f) + end
        out.append(line)
    dst.write_text("".join(out))


def make_conf(run):
    src = CONF_DIR / f"analyzer_{run_tag(run)}_d5_{SUFFIX}.conf"
    dst = CONF_DIR / f"analyzer_{run_tag(run)}_d5_{SUFFIX}_{TAG}.conf"
    text = src.read_text()
    geo = f"param/DCGEO/{config.SUB_DIR}/DCGeomParam_{run_tag(run)}_{SUFFIX}"
    if geo + "\n" not in text:
        raise RuntimeError(f"DCGEO line not found in {src}")
    dst.write_text(text.replace(geo + "\n", f"{geo}_{TAG}\n"))
    return dst


def make_runlists(runs):
    head, body, body_tag = None, [], []
    for run in runs:
        src = RUNLIST_DIR / config.SUB_DIR / f"d5_{run_tag(run)}.yml"
        text = src.read_text()
        h, b = text.split("RUN:\n", 1)
        head = head or h + "RUN:\n"
        body.append(b)
        b_tag = (b.replace(f"_{SUFFIX}_d5:", f"_{SUFFIX}_d5{TAG}:")
                  .replace(f"_d5_{SUFFIX}.conf", f"_d5_{SUFFIX}_{TAG}.conf")
                  .replace(f"_D5_{SUFFIX}.root", f"_D5_{SUFFIX}_{TAG}.root"))
        body_tag.append(b_tag)
    label = runs_label(runs)
    p_norm = RUNLIST_DIR / "scratch" / f"d5_{label}.yml"
    p_tag = RUNLIST_DIR / "scratch" / f"d5{TAG}_{label}.yml"
    p_norm.write_text(head + "".join(body))
    p_tag.write_text(head + "".join(body_tag))
    return p_norm, p_tag


def do_prep(args):
    for run in args.runs:
        inflate_blc1_res(dcgeo_path(run), Path(f"{dcgeo_path(run)}_{TAG}"), args.factor)
        make_conf(run)
        print(colored(f"[INFO] {run_tag(run)}: DCGEO_{TAG} (BLC1 Res x{args.factor:g}) and conf written", "green"))
    p_norm, p_tag = make_runlists(args.runs)
    print(colored("[INFO] move existing D5 outputs aside if they should be kept, then submit:", "yellow"))
    for p in (p_norm, p_tag):
        print(f"  cd {config.ANALYZER_DIR}/runmanager && ./myrun.py runlist/scratch/{p.name} --local")


# -- review -------------------------------------------------------------------
def parse_kv(path):
    kv = {}
    for line in path.read_text().splitlines():
        p = line.split(None, 1)
        if len(p) == 2:
            kv[p[0]] = p[1]
    return kv


def review_one(run):
    res = {"run": run}
    rc, _ = run_cmd(f"{BIN_DIR}/d5_bc1reso_review {run} {BEAM_FLAG} --d5 {d5_root(run, TAG)} "
                    f"--dcgeo {dcgeo_path(run)}", capture=True)
    txt_y = txt_dir(run) / f"d5_bc1reso_review_{run_tag(run)}.txt"
    if rc == 0 and txt_y.exists():
        kv = parse_kv(txt_y)
        res["dy_a"], res["dy_b"] = float(kv["dy_a_mm"]), float(kv["dy_b_mm"])
        res["dy_err"] = max(float(kv["dy_a_err_mm"]), float(kv["dy_b_err_mm"]))
    rc2, _ = run_cmd(f"{BIN_DIR}/d5_fitfree_review {run} {BEAM_FLAG} --d5 {d5_root(run)} "
                     f"--dcgeo {dcgeo_path(run)}", capture=True)
    txt_x = txt_dir(run) / f"d5_fitfree_review_{run_tag(run)}.txt"
    if rc2 == 0 and txt_x.exists():
        kv = parse_kv(txt_x)
        res["s"] = float(kv["blc1_ab_opposite_shift_s_mm"])
        res["ddxu"] = float(kv["ddxu_mean_percent"])
        res["dv"] = float(kv["dv_mean_mrad"])
    return res


def fmt(v, f="{:+.4f}"):
    return "n/a" if v is None else f.format(v)


def do_review(args):
    rows = [review_one(run) for run in args.runs]
    lines = ["run      dy_a[mm]  dy_b[mm]  err[mm]   s[mm]     ddxu[%]   dv[mrad]  converged"]
    for r in rows:
        ok = all(v is not None for v in (r.get("dy_a"), r.get("dy_b"), r.get("s")))
        tol_y = max(args.tol, 2.0 * r.get("dy_err", 0.0))
        conv = ok and abs(r["dy_a"]) < tol_y and abs(r["dy_b"]) < tol_y and abs(r["s"]) < args.tol
        lines.append(f"{r['run']:05d}    {fmt(r.get('dy_a'))}   {fmt(r.get('dy_b'))}   {fmt(r.get('dy_err'))}   {fmt(r.get('s'))}   "
                     f"{fmt(r.get('ddxu'))}   {fmt(r.get('dv'))}   {'yes' if conv else 'no'}")
    out_dir = SUMMARY_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = out_dir / f"d5_align_summary_{runs_label(args.runs)}.txt"
    summary.write_text("\n".join(lines) + f"\n(tolerance {args.tol} mm, dy: max(tol, 2 err); dy from d5_bc1reso_review {BEAM_FLAG}, s from d5_fitfree_review {BEAM_FLAG})\n")
    print("\n".join(lines))
    print(colored(f"[INFO] summary: {summary}", "green"))


# -- apply --------------------------------------------------------------------
def do_apply(args):
    use_y = args.y or not args.x
    use_x = args.x or not args.y
    backup = None
    if args.write:
        backup = BACKUP_BASE / f"dcgeo_d5align_{datetime.datetime.now():%Y%m%d_%H%M%S}"
        backup.mkdir(parents=True, exist_ok=False)
    for run in args.runs:
        dx_a = dx_b = dy_a = dy_b = 0.0
        if use_y:
            kv = parse_kv(txt_dir(run) / f"d5_bc1reso_review_{run_tag(run)}.txt")
            dy_a, dy_b = float(kv["dy_a_mm"]), float(kv["dy_b_mm"])
        if use_x:
            s = float(parse_kv(txt_dir(run) / f"d5_fitfree_review_{run_tag(run)}.txt")["blc1_ab_opposite_shift_s_mm"])
            dx_a, dx_b = s, -s
        opt = (f"--skip-xy --residual-offset BLC1a:{dx_a:.5f},{dy_a:.5f} "
               f"--residual-offset BLC1b:{dx_b:.5f},{dy_b:.5f}")
        if backup:
            shutil.copy2(dcgeo_path(run), backup / dcgeo_path(run).name)
            opt += " --write"
        rc, out = run_cmd(f"python3 {APPLY_SCRIPT} {dcgeo_path(run)} {opt}", capture=True)
        print(out.strip().splitlines()[-1] if out.strip() else "")
        if rc != 0:
            print(colored(f"[Error] apply failed for {run_tag(run)}", "red"))
            sys.exit(rc)
    if backup:
        print(colored(f"[INFO] DCGEO backup: {backup}", "green"))
        dcprm_flag = " --kaon" if SUFFIX == "K" else ""
        print(colored(f"[INFO] next: myrun BcIn -> run_dcprm.py <run> resi --bcin{dcprm_flag} -> run_d5align.py prep", "yellow"))
    else:
        print(colored("[INFO] dry-run only (add --write)", "yellow"))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mode", choices=["prep", "review", "apply"])
    parser.add_argument("runs", type=int, nargs="+", help="run numbers")
    parser.add_argument("--factor", type=float, default=1000.0, help="prep: BLC1 Res multiplier (default 1000)")
    parser.add_argument("--tol", type=float, default=0.01, help="review: convergence tolerance [mm] (default 0.01)")
    parser.add_argument("--y", action="store_true", help="apply: y (BLC2-referenced dy) only")
    parser.add_argument("--x", action="store_true", help="apply: x (BLC1a/b opposite shift) only")
    parser.add_argument("--write", action="store_true", help="apply: back up and write DCGEO")
    parser.add_argument("--beam", choices=sorted(BEAM), default="kaon",
                        help="beam used for alignment: review selects it, params use its suffix (default kaon)")
    args = parser.parse_args()
    global SUFFIX, BEAM_FLAG
    SUFFIX, BEAM_FLAG = BEAM[args.beam]
    {"prep": do_prep, "review": do_review, "apply": do_apply}[args.mode](args)


if __name__ == "__main__":
    main()
