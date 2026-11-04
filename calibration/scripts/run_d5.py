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
  pos  Measure BLC1a/b dy (BLC2-referenced fit) and the a/b opposite x shift s
       (<delta_x - delta_u>) with d5_pos_review and print the table.
       Not converged: write them into the BLC1a/b Ofs of the DCGEO (backup first).
       Converged: write nothing and make the check PDF (d5_plane_resi_qa).
       --debug: measure and print only.

One cycle:
  ref -> myrun D5 (normal, BLC2-referenced) -> pos -> myrun BcIn -> run_dcprm.py <run> resi --bcin
  -> ref -> ... until pos reports converged.
Score: max(|dy_a|, |dy_b|) / max(--tol, 2 * statistical error) and |s| / --tol, the larger one.
Converged: score < 1.
Settled: outside the criterion, and the score (below) did not shrink below POS_IMPROVE of
the previous pass in the last 2 of 3 passes (history OUTPUT_DIR/calib/runXXXXX/d5_pos_<S>.txt,
calib_convergence.py): nothing written, a person looks.
Exit code of pos (lib/runner.py, worst run): 3 some run settled, 2 Ofs written for some run,
0 all converged (or --debug), 1 failed.
The common x shift of BLC1 is degenerate with delta and is left unchanged.
"""

import re
import sys
from pathlib import Path
from termcolor import colored

sys.path.append(str(Path(__file__).resolve().parents[2]))
from lib import config, runner  # noqa: E402
import calib_convergence as conv  # noqa: E402

DCGEO_DIR = config.PARAM_DIR / "DCGEO" / config.SUB_DIR
CONF_DIR = config.PARAM_DIR / "conf" / config.SUB_DIR
RUNLIST_DIR = config.ANALYZER_DIR / "runmanager" / "runlist"
BIN_DIR = runner.BIN_DIR
APPLY_SCRIPT = runner.SCRIPT_DIR / "apply_dcgeo_offset.py"
SUMMARY_DIR = config.WORK_DIR / "d5_align"  # number files (not in results/img)
REF_DIR = SUMMARY_DIR / "ref"  # throw-away DCGEO copy and conf of the BLC2-referenced fit
CONF_FILE_KEYS = ("UNPACK", "DIGIT", "CMAP", "USER", "HDPRM", "HDPHC", "DCGEO", "DCTDC", "DCDRFT",
                  "TPCPRM", "TPCPOS", "FLDMAP", "D5MTX")
TAG = "bc1reso"
POS_IMPROVE = 0.8  # a pos pass improves when its score shrinks below this fraction of the previous one


def run_tag(run):
    return f"run{run:05d}"


def beam_flag(s):
    return "--kaon" if s == "K" else "--pion"


def kaon_opt(s):
    return " --kaon" if s == "K" else ""


def dcgeo_path(run, s):
    return DCGEO_DIR / f"DCGeomParam_{run_tag(run)}_{s}"


def ref_dcgeo_path(run, s):
    return REF_DIR / run_tag(run) / f"DCGeomParam_{run_tag(run)}_{s}_{TAG}"


def ref_conf_path(run, s):
    return REF_DIR / run_tag(run) / f"analyzer_{run_tag(run)}_d5_{s}_{TAG}.conf"


def d5_root(run, s, tag=None):
    """D5 normal in DECODE_DIR; the BLC2-referenced (throw-away) one in SCRATCH_DIR."""
    if tag:
        return config.SCRATCH_DIR / run_tag(run) / f"{run_tag(run)}_D5_{s}_{tag}.root"
    return config.DECODE_DIR / run_tag(run) / f"{run_tag(run)}_D5_{s}.root"


def chamber_root(run, det, s):
    return config.DECODE_DIR / run_tag(run) / f"{run_tag(run)}_{det}_resi_{s}.root"


def txt_dir(run):
    return SUMMARY_DIR / run_tag(run)


def runs_label(runs, s):
    label = f"{runs[0]:05d}-{runs[-1]:05d}_{len(runs)}runs" if len(runs) > 1 else f"{runs[0]:05d}"
    return label + ("" if s == "K" else f"_{s}")


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


def make_conf(run, s):
    """Copy of the D5 conf with DCGEO -> the Res-inflated copy; every file path made absolute."""
    src = CONF_DIR / f"analyzer_{run_tag(run)}_d5_{s}.conf"
    if not src.exists():
        runner.fail(f"{src} not found (make the D5 conf of suffix {s} first)")
    out, found = [], False
    for line in src.read_text().splitlines():
        f = line.split()
        key = f[0].rstrip(":") if f else ""
        if key in CONF_FILE_KEYS and len(f) >= 2:
            path = Path(f[1]) if Path(f[1]).is_absolute() else config.ANALYZER_DIR / f[1]
            if key == "DCGEO":
                if path.resolve() != dcgeo_path(run, s).resolve():
                    runner.fail(f"{src}: DCGEO is {f[1]}, expected {dcgeo_path(run, s)}")
                path, found = ref_dcgeo_path(run, s), True
            line = f"{key}:\t{path}"
        out.append(line)
    if not found:
        runner.fail(f"no DCGEO line in {src}")
    ref_conf_path(run, s).write_text("\n".join(out) + "\n")


def ref_runlist_entry(run, s, body):
    """The run's D5 runlist entry pointed at the Res-inflated conf and the SCRATCH_DIR output."""
    conf_rel = f"param/conf/{config.SUB_DIR}/analyzer_{run_tag(run)}_d5_{s}.conf"
    root_abs = str(d5_root(run, s))
    if f"_{s}_d5:" not in body or conf_rel not in body or root_abs not in body:
        runner.fail(f"d5_{run_tag(run)}.yml: no {s} D5 entry with conf {conf_rel} and root {root_abs}")
    return (body.replace(f"_{s}_d5:", f"_{s}_d5{TAG}:")
            .replace(conf_rel, str(ref_conf_path(run, s)))
            .replace(root_abs, str(d5_root(run, s, TAG))))


def make_runlists(runs, s):
    head, body, body_tag = None, [], []
    for run in runs:
        h, b = (RUNLIST_DIR / config.SUB_DIR / f"d5_{run_tag(run)}.yml").read_text().split("RUN:\n", 1)
        head = head or h + "RUN:\n"
        body.append(b)
        body_tag.append(ref_runlist_entry(run, s, b))
    out_dir = RUNLIST_DIR / "scratch" / "d5"
    out_dir.mkdir(parents=True, exist_ok=True)
    p_norm = out_dir / f"d5_{runs_label(runs, s)}.yml"
    p_tag = out_dir / f"d5{TAG}_{runs_label(runs, s)}.yml"
    p_norm.write_text(head + "".join(body))
    p_tag.write_text(head + "".join(body_tag))
    return p_norm, p_tag


def do_ref(args, s):
    for run in args.runs:
        if not dcgeo_path(run, s).exists():
            runner.fail(f"{dcgeo_path(run, s)} not found")
        ref_dcgeo_path(run, s).parent.mkdir(parents=True, exist_ok=True)
        d5_root(run, s, TAG).parent.mkdir(parents=True, exist_ok=True)
        inflate_blc1_res(dcgeo_path(run, s), ref_dcgeo_path(run, s), args.factor)
        make_conf(run, s)
        runner.info(f"{run_tag(run)}: {ref_dcgeo_path(run, s)} (BLC1 Res x{args.factor:g}) and its conf written")
    p_norm, p_tag = make_runlists(args.runs, s)
    runner.warn("decode both (move existing D5 outputs aside first if they should be kept):")
    for p in (p_norm, p_tag):
        print(f"  cd {config.ANALYZER_DIR}/runmanager && ./myrun.py {p.relative_to(RUNLIST_DIR.parent)} --local")
    runner.warn(f"then: run_d5.py {' '.join(map(str, args.runs))} pos{kaon_opt(s)}")


# -- pos ----------------------------------------------------------------------
def copy_matches(run, s):
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
    copy = ref_dcgeo_path(run, s)
    return copy.exists() and rows(copy) == rows(dcgeo_path(run, s))


def parse_kv(path):
    kv = {}
    for line in path.read_text().splitlines():
        p = line.split(None, 1)
        if len(p) == 2:
            kv[p[0]] = p[1]
    return kv


def measure(run, s):
    """dy_a, dy_b, dy_err, s, ddxu, dv of d5_pos_review; None if it failed."""
    rc, _ = runner.sh(f"{BIN_DIR}/d5_pos_review {run} {beam_flag(s)} --d5 {d5_root(run, s)} "
                      f"--d5ref {d5_root(run, s, TAG)} --dcgeo {dcgeo_path(run, s)}", capture=True, check=False)
    txt = txt_dir(run) / f"d5_pos_review_{run_tag(run)}_{s}.txt"
    if rc != 0 or not txt.exists():
        return None
    kv = parse_kv(txt)
    keys = {"dy_a": "dy_a_mm", "dy_b": "dy_b_mm", "dy_err_a": "dy_a_err_mm", "dy_err_b": "dy_b_err_mm",
            "s": "blc1_ab_opposite_shift_s_mm", "ddxu": "ddxu_mean_percent", "dv": "dv_mean_mrad"}
    try:
        r = {k: float(kv[v]) for k, v in keys.items()}
    except (KeyError, ValueError):
        runner.warn(f"{txt}: missing or bad values")
        return None
    r["dy_err"] = max(r.pop("dy_err_a"), r.pop("dy_err_b"))
    return r


def pos_score(r, tol):
    """Largest of |dy_a|, |dy_b|, |s| in units of their tolerances (< 1 when converged)."""
    tol_y = max(tol, 2.0 * r["dy_err"])
    return max(abs(r["dy_a"]) / tol_y, abs(r["dy_b"]) / tol_y, abs(r["s"]) / tol)


def judge_pos(run, r, s, tol):
    path = config.OUTPUT_DIR / "calib" / run_tag(run) / f"d5_pos_{s}.txt"
    details = f"dy_a={r['dy_a']:+.5f} dy_b={r['dy_b']:+.5f} err={r['dy_err']:.5f} s={r['s']:+.5f}"
    score = pos_score(r, tol)
    return conv.judge_at(path, "BLC1", score, score < 1.0, lambda prev, new: new > POS_IMPROVE * prev, details)


def fmt(v, f="{:+.4f}"):
    return "n/a" if v is None else f.format(v)


def summary_table(rows, s, tol):
    """rows: {run: measurement or None}."""
    lines = ["run      dy_a[mm]  dy_b[mm]  err[mm]   s[mm]     ddxu[%]   dv[mrad]  state"]
    for run, r in rows.items():
        r = r or {}
        state = r.get("state") or ("n/a" if not r else "yes" if pos_score(r, tol) < 1.0 else "no")
        lines.append(f"{run:05d}    {fmt(r.get('dy_a'))}   {fmt(r.get('dy_b'))}   {fmt(r.get('dy_err'))}   "
                     f"{fmt(r.get('s'))}   {fmt(r.get('ddxu'))}   {fmt(r.get('dv'))}   {state}")
    return "\n".join(lines) + f"\n(beam {beam_flag(s)}; converged: |dy| < max({tol} mm, 2 err), |s| < {tol} mm)\n"


def write_summary(rows, s, tol):
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    path = SUMMARY_DIR / f"d5_pos_summary_{runs_label(list(rows), s)}.txt"
    path.write_text(summary_table(rows, s, tol))
    print(summary_table(rows, s, tol), end="")
    runner.info(f"summary: {path}")


def write_ofs(run, r, s):
    runner.backup_file(dcgeo_path(run, s))
    opt = (f"--skip-xy --residual-offset BLC1a:{r['s']:.5f},{r['dy_a']:.5f} "
           f"--residual-offset BLC1b:{-r['s']:.5f},{r['dy_b']:.5f} --write")
    _, out = runner.sh(f"python3 {APPLY_SCRIPT} {dcgeo_path(run, s)} {opt}", capture=True)
    print(out.strip().splitlines()[-1] if out.strip() else "")


def make_check_pdf(run, s):
    runner.sh(f"{BIN_DIR}/d5_plane_resi_qa {run} {beam_flag(s)} --d5 {d5_root(run, s)} "
              f"--bcin {chamber_root(run, 'BcIn', s)} --bcout {chamber_root(run, 'BcOut', s)}", check=False)


def check_pos_inputs(run, s, force):
    if not dcgeo_path(run, s).exists():
        runner.fail(f"{dcgeo_path(run, s)} not found (suffix {s}; --kaon for K)")
    if not copy_matches(run, s):
        runner.fail(f"{ref_dcgeo_path(run, s)} does not match {dcgeo_path(run, s).name} "
                    f"(DCGEO changed after ref): run 'run_d5.py {run} ref' and decode again")
    for root in (d5_root(run, s), d5_root(run, s, TAG)):
        runner.checked_root(root, s, force)


def do_pos(args, s):
    for run in args.runs:
        check_pos_inputs(run, s, args.force)
    rows = {run: measure(run, s) for run in args.runs}
    failed = [run for run, r in rows.items() if r is None]
    if not (failed or args.debug):
        for run, r in rows.items():
            r["state"] = judge_pos(run, r, s, args.tol)
    write_summary(rows, s, args.tol)
    if failed:
        runner.fail(f"measurement failed for {failed} (n/a); nothing written")
    if args.debug:
        runner.warn("--debug: nothing written")
        return

    written = [run for run, r in rows.items() if r["state"] == conv.WRITTEN]
    for run in written:
        write_ofs(run, rows[run], s)
    if written:
        rl = " ".join(map(str, written))
        runner.warn(f"next for {rl}: myrun BcIn -> run_dcprm.py <run> resi --bcin{kaon_opt(s)} -> run_d5.py {rl} ref")
    for run, r in rows.items():
        if r["state"] == conv.SETTLED:
            print(colored(f"[SETTLED] {run_tag(run)}: not converging over {conv.SETTLE_PASSES} passes; "
                          "nothing written, look at the d5_pos_review PDF", "red", attrs=["bold"]))
        elif r["state"] == conv.CONVERGED:
            runner.info(f"{run_tag(run)} converged: check PDF")
            make_check_pdf(run, s)
    sys.exit(conv.exit_code([r["state"] for r in rows.values()]))


def main():
    parser = runner.base_parser(__doc__, ["ref", "pos"])
    parser.add_argument("--factor", type=float, default=1000.0, help="ref: BLC1 Res multiplier (default 1000)")
    parser.add_argument("--tol", type=float, default=0.01, help="pos: convergence tolerance [mm] (default 0.01)")
    args = parser.parse_args()
    {"ref": do_ref, "pos": do_pos}[args.mode](args, runner.suffix_of(args))


if __name__ == "__main__":
    main()
