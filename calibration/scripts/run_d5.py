#!/usr/bin/env python3
"""BLC1 position against BLC2 through the D5 matrix (BLC2 fixed, BLC1a and BLC1b translated).

Usage: run_d5.py <run...> {ref,pos} [--kaon] [--debug] [--force]
  Several runs are processed one by one. Particle: pi by default, K with --kaon (param suffix
  _Pi / _K, beam flag of the D5 tools). This script prepares params and fits only; the decodes
  (myrun) are run outside.

  ref  Make the inputs of the BLC2-referenced D5 fit: a copy of the DCGEO with the BLC1a/b Res
       multiplied by --factor (the fit is then set by BLC2) and its conf (absolute paths), both
       throw-away files in REF_DIR/runXXXXX/, and the runlists in runmanager/runlist/scratch/d5/,
       and print the myrun commands (D5 normal; D5 BLC2-referenced, written to SCRATCH_DIR).
  pos  Measure BLC1a/b dy (d5_bc1reso_review, BLC2-referenced fit) and the a/b opposite x shift
       s (d5_fitfree_review, <delta_x - delta_u>) and print the table.
       Not converged: write them into the BLC1a/b Ofs of the DCGEO (backup first).
       Converged: write nothing and make the check PDF (d5_plane_resi_qa).
       --debug: measure and print only.

One cycle:
  ref -> myrun D5 (normal, BLC2-referenced) -> pos -> myrun BcIn -> run_dcprm.py <run> resi --bcin
  -> ref -> ... until pos reports converged.
Converged: |dy| < max(--tol, 2 * statistical error) and |s| < --tol.
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
from lib.input_check import check_input

DCGEO_DIR = config.PARAM_DIR / "DCGEO" / config.SUB_DIR
CONF_DIR = config.PARAM_DIR / "conf" / config.SUB_DIR
RUNLIST_DIR = config.ANALYZER_DIR / "runmanager" / "runlist"
BIN_DIR = project_root / "bin"
APPLY_SCRIPT = Path(__file__).resolve().parent / "apply_dcgeo_offset.py"
BACKUP_BASE = Path("/group/had/sks/Users/sryuta/backup")
SUMMARY_DIR = Path("/group/had/sks/Users/sryuta/tmp/d5_align")  # number files (not in results/img)
REF_DIR = SUMMARY_DIR / "ref"  # throw-away DCGEO copy and conf of the BLC2-referenced fit
CONF_FILE_KEYS = ("UNPACK", "DIGIT", "CMAP", "USER", "HDPRM", "HDPHC", "DCGEO", "DCTDC", "DCDRFT",
                  "TPCPRM", "TPCPOS", "FLDMAP", "D5MTX")
TAG = "bc1reso"

# set from --kaon in main()
SUFFIX = "Pi"
BEAM_FLAG = "--pion"


def run_tag(run):
    return f"run{run:05d}"


def dcgeo_path(run):
    return DCGEO_DIR / f"DCGeomParam_{run_tag(run)}_{SUFFIX}"


def ref_dcgeo_path(run):
    return REF_DIR / run_tag(run) / f"DCGeomParam_{run_tag(run)}_{SUFFIX}_{TAG}"


def ref_conf_path(run):
    return REF_DIR / run_tag(run) / f"analyzer_{run_tag(run)}_d5_{SUFFIX}_{TAG}.conf"


def d5_root(run, tag=None):
    """D5 normal in DECODE_DIR; the BLC2-referenced (throw-away) one in SCRATCH_DIR."""
    if tag:
        return config.SCRATCH_DIR / run_tag(run) / f"{run_tag(run)}_D5_{SUFFIX}_{tag}.root"
    return config.DECODE_DIR / run_tag(run) / f"{run_tag(run)}_D5_{SUFFIX}.root"


def chamber_root(run, det):
    return config.DECODE_DIR / run_tag(run) / f"{run_tag(run)}_{det}_resi_{SUFFIX}.root"


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


def is_blc1_row(fields):
    return len(fields) >= 13 and fields[0].isdigit() and re.match(r"BLC1[ab]-", fields[1])


# -- ref ----------------------------------------------------------------------
def inflate_blc1_res(src, dst, factor):
    out = []
    for line in src.read_text().splitlines(keepends=True):
        f = line.split()
        if is_blc1_row(f):
            f[9] = f"{float(f[9]) * factor:.6g}"
            end = "\n" if line.endswith("\n") else ""
            line = "\t".join(f) + end
        out.append(line)
    dst.write_text("".join(out))


def make_conf(run):
    """Copy of the D5 conf with DCGEO -> the Res-inflated copy; every file path made absolute."""
    src = CONF_DIR / f"analyzer_{run_tag(run)}_d5_{SUFFIX}.conf"
    if not src.exists():
        raise RuntimeError(f"{src} not found (make the D5 conf of suffix {SUFFIX} first)")
    out, found = [], False
    for line in src.read_text().splitlines():
        f = line.split()
        key = f[0].rstrip(":") if f else ""
        if key in CONF_FILE_KEYS and len(f) >= 2:
            path = Path(f[1]) if Path(f[1]).is_absolute() else config.ANALYZER_DIR / f[1]
            if key == "DCGEO":
                if path.resolve() != dcgeo_path(run).resolve():
                    raise RuntimeError(f"{src}: DCGEO is {f[1]}, expected {dcgeo_path(run)}")
                path, found = ref_dcgeo_path(run), True
            line = f"{key}:\t{path}"
        out.append(line)
    if not found:
        raise RuntimeError(f"no DCGEO line in {src}")
    dst = ref_conf_path(run)
    dst.write_text("\n".join(out) + "\n")
    return dst


def make_runlists(runs):
    head, body, body_tag = None, [], []
    for run in runs:
        src = RUNLIST_DIR / config.SUB_DIR / f"d5_{run_tag(run)}.yml"
        text = src.read_text()
        h, b = text.split("RUN:\n", 1)
        if f"_{SUFFIX}_d5:" not in b or f"_D5_{SUFFIX}.root" not in b:
            raise RuntimeError(f"{src} has no {SUFFIX} D5 entry")
        head = head or h + "RUN:\n"
        body.append(b)
        conf_rel = f"param/conf/{config.SUB_DIR}/analyzer_{run_tag(run)}_d5_{SUFFIX}.conf"
        root_abs = str(d5_root(run))
        if conf_rel not in b or root_abs not in b:
            raise RuntimeError(f"{src}: conf {conf_rel} / root {root_abs} not found")
        body_tag.append(b.replace(f"_{SUFFIX}_d5:", f"_{SUFFIX}_d5{TAG}:")
                         .replace(conf_rel, str(ref_conf_path(run)))
                         .replace(root_abs, str(d5_root(run, TAG))))
    label = runs_label(runs) + ("" if SUFFIX == "K" else f"_{SUFFIX}")
    out_dir = RUNLIST_DIR / "scratch" / "d5"
    out_dir.mkdir(parents=True, exist_ok=True)
    p_norm = out_dir / f"d5_{label}.yml"
    p_tag = out_dir / f"d5{TAG}_{label}.yml"
    p_norm.write_text(head + "".join(body))
    p_tag.write_text(head + "".join(body_tag))
    return p_norm, p_tag


def do_ref(args):
    for run in args.runs:
        if not dcgeo_path(run).exists():
            print(colored(f"[Error] {dcgeo_path(run)} not found", "red"))
            sys.exit(1)
        ref_dcgeo_path(run).parent.mkdir(parents=True, exist_ok=True)
        d5_root(run, TAG).parent.mkdir(parents=True, exist_ok=True)
        inflate_blc1_res(dcgeo_path(run), ref_dcgeo_path(run), args.factor)
        make_conf(run)
        print(colored(f"[INFO] {run_tag(run)}: {ref_dcgeo_path(run)} (BLC1 Res x{args.factor:g}) "
                      "and its conf written", "green"))
    p_norm, p_tag = make_runlists(args.runs)
    print(colored("[INFO] decode both (move existing D5 outputs aside first if they should be kept):", "yellow"))
    for p in (p_norm, p_tag):
        print(f"  cd {config.ANALYZER_DIR}/runmanager && ./myrun.py {p.relative_to(RUNLIST_DIR.parent)} --local")
    print(colored(f"[INFO] then: run_d5.py {' '.join(map(str, args.runs))} pos"
                  + (" --kaon" if SUFFIX == "K" else ""), "yellow"))


# -- pos ----------------------------------------------------------------------
def copy_matches(run):
    """The Res-inflated copy differs from the DCGEO only in the BLC1 Res column."""
    def rows(path):
        out = {}
        for line in path.read_text().splitlines():
            f = line.split()
            if len(f) >= 13 and f[0].isdigit():
                if is_blc1_row(f):
                    f = f[:9] + f[10:]
                out[f[0]] = f
        return out
    copy = ref_dcgeo_path(run)
    return copy.exists() and rows(copy) == rows(dcgeo_path(run))


def parse_kv(path):
    kv = {}
    for line in path.read_text().splitlines():
        p = line.split(None, 1)
        if len(p) == 2:
            kv[p[0]] = p[1]
    return kv


def measure(run):
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


def is_converged(r, tol):
    if any(r.get(k) is None for k in ("dy_a", "dy_b", "s")):
        return None
    tol_y = max(tol, 2.0 * r.get("dy_err", 0.0))
    return abs(r["dy_a"]) < tol_y and abs(r["dy_b"]) < tol_y and abs(r["s"]) < tol


def fmt(v, f="{:+.4f}"):
    return "n/a" if v is None else f.format(v)


def write_ofs(rows):
    backup = BACKUP_BASE / f"dcgeo_d5pos_{datetime.datetime.now():%Y%m%d_%H%M%S}"
    backup.mkdir(parents=True, exist_ok=False)
    for r in rows:
        run = r["run"]
        shutil.copy2(dcgeo_path(run), backup / dcgeo_path(run).name)
        opt = (f"--skip-xy --residual-offset BLC1a:{r['s']:.5f},{r['dy_a']:.5f} "
               f"--residual-offset BLC1b:{-r['s']:.5f},{r['dy_b']:.5f} --write")
        rc, out = run_cmd(f"python3 {APPLY_SCRIPT} {dcgeo_path(run)} {opt}", capture=True)
        print(out.strip().splitlines()[-1] if out.strip() else "")
        if rc != 0:
            print(colored(f"[Error] writing Ofs failed for {run_tag(run)}", "red"))
            sys.exit(rc)
    print(colored(f"[INFO] DCGEO backup: {backup}", "green"))


def make_check_pdf(run):
    flag = "--kaon" if SUFFIX == "K" else "--pion"
    run_cmd(f"{BIN_DIR}/d5_plane_resi_qa {run} {flag} --d5 {d5_root(run)} "
            f"--bcin {chamber_root(run, 'BcIn')} --bcout {chamber_root(run, 'BcOut')}")


def do_pos(args):
    for run in args.runs:
        if not dcgeo_path(run).exists():
            print(colored(f"[Error] {dcgeo_path(run)} not found (suffix {SUFFIX}; --kaon for K)", "red"))
            sys.exit(1)
        if not copy_matches(run):
            print(colored(f"[Error] {ref_dcgeo_path(run)} does not match {dcgeo_path(run).name} "
                          f"(DCGEO changed after ref): run 'run_d5.py {run} ref' and decode again", "red"))
            sys.exit(1)
        for root in (d5_root(run), d5_root(run, TAG)):
            if not root.exists():
                print(colored(f"[Error] {root} not found", "red"))
                sys.exit(1)
            if not check_input(root, SUFFIX, args.force):
                sys.exit(1)

    rows = [measure(run) for run in args.runs]
    lines = ["run      dy_a[mm]  dy_b[mm]  err[mm]   s[mm]     ddxu[%]   dv[mrad]  converged"]
    for r in rows:
        conv = is_converged(r, args.tol)
        r["converged"] = conv
        lines.append(f"{r['run']:05d}    {fmt(r.get('dy_a'))}   {fmt(r.get('dy_b'))}   {fmt(r.get('dy_err'))}   "
                     f"{fmt(r.get('s'))}   {fmt(r.get('ddxu'))}   {fmt(r.get('dv'))}   "
                     f"{'n/a' if conv is None else 'yes' if conv else 'no'}")
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    label = runs_label(args.runs) + ("" if SUFFIX == "K" else f"_{SUFFIX}")
    summary = SUMMARY_DIR / f"d5_pos_summary_{label}.txt"
    summary.write_text("\n".join(lines) + f"\n(beam {BEAM_FLAG}; converged: |dy| < max({args.tol} mm, 2 err), "
                       f"|s| < {args.tol} mm)\n")
    print("\n".join(lines))
    print(colored(f"[INFO] summary: {summary}", "green"))

    if any(r["converged"] is None for r in rows):
        print(colored("[Error] measurement failed for some runs (n/a); nothing written", "red"))
        sys.exit(1)
    todo = [r for r in rows if not r["converged"]]
    done = [r for r in rows if r["converged"]]
    if args.debug:
        print(colored("[INFO] --debug: nothing written", "yellow"))
        return
    if todo:
        write_ofs(todo)
        rl = " ".join(str(r["run"]) for r in todo)
        print(colored(f"[INFO] next for {rl}: myrun BcIn -> run_dcprm.py <run> resi --bcin"
                      + (" --kaon" if SUFFIX == "K" else "") + f" -> run_d5.py {rl} ref", "yellow"))
    for r in done:
        print(colored(f"[INFO] {run_tag(r['run'])} converged: check PDF", "green"))
        make_check_pdf(r["run"])


def main():
    global SUFFIX, BEAM_FLAG
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", type=int, nargs="+", help="run numbers (processed one by one)")
    parser.add_argument("mode", choices=["ref", "pos"])
    parser.add_argument("--kaon", action="store_true", help="K (param suffix _K, K window); default pi")
    parser.add_argument("--factor", type=float, default=1000.0, help="ref: BLC1 Res multiplier (default 1000)")
    parser.add_argument("--tol", type=float, default=0.01, help="pos: convergence tolerance [mm] (default 0.01)")
    parser.add_argument("--debug", action="store_true", help="pos: measure and print only")
    parser.add_argument("--force", action="store_true",
                        help="pos: run even if a D5 decode is older than its params")
    args = parser.parse_args()
    if args.kaon:
        SUFFIX, BEAM_FLAG = "K", "--kaon"
    {"ref": do_ref, "pos": do_pos}[args.mode](args)


if __name__ == "__main__":
    main()
