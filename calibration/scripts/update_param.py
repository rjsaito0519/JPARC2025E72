#!/usr/bin/env python3
"""Write the fit results of run_hodo.py / run_dcprm.py into the per-run param files.

Usage: update_param.py <run> {K,Pi} <ptype> [--bcin|--bcout|--all] [--separate] [--ftof]
  hdprm, t0, hdphc   hodoscopes (detector set of lib/config.py HODO_DETECTORS; --ftof the FTOF set)
  htofprm, htofphc   HTOF_Calib -> HDPRM / HDPHC
  dctdc              BLC_TDC -> DCTDC
  residual           BLC_residual -> DCGEO Ofs (the correction is added to the current Ofs)
  resolution         BLC_pull -> DCGEO Res: Res_new = p * mean(Res in group) for every layer of the group
                     (default one group per chamber, --separate a and b, --all BLC1+BLC2 with p averaged)
Chamber: --bcout (default) or --bcin.
Fit outputs: OUTPUT_DIR/root/runXXXXX/runXXXXX_<det>[_<kind>]_<S>.root. A fit output older than its decode
stops everything (exit 1). A missing hodo output is skipped (its rows keep their values); a missing DC
output stops (exit 1).
residual and resolution judge each group (calib_convergence.py) and write only the groups still moving;
exit code 0 converged, 2 written, 3 settled (lib/runner.py).
"""
import re
import sys
from pathlib import Path

import numpy as np
from termcolor import colored

sys.path.append(str(Path(__file__).resolve().parents[2]))
from lib import config, runner  # noqa: E402
from lib.param import ParamFile  # noqa: E402
import calib_convergence as conv  # noqa: E402
import phc_conf  # noqa: E402
import update_dctdc  # noqa: E402
import update_hdprm  # noqa: E402
import update_phc  # noqa: E402
import update_residual  # noqa: E402

CATEGORIES = {
    "hdprm":      {"dir": "HDPRM", "tpl": "HodoParam_e72_example",    "prefix": "HodoParam",    "key_len": 5, "start_col": 5},
    "htofprm":    {"dir": "HDPRM", "tpl": "HodoParam_e72_example",    "prefix": "HodoParam",    "key_len": 5, "start_col": 5},
    "t0":         {"dir": "HDPRM", "tpl": "HodoParam_e72_example",    "prefix": "HodoParam",    "key_len": 5, "start_col": 5},
    "hdphc":      {"dir": "HDPHC", "tpl": "HodoPHCParam_e72_example", "prefix": "HodoPHCParam", "key_len": 4, "start_col": 4},
    "htofphc":    {"dir": "HDPHC", "tpl": "HodoPHCParam_e72_example", "prefix": "HodoPHCParam", "key_len": 4, "start_col": 4},
    "dctdc":      {"dir": "DCTDC", "tpl": "DCTdcParam_e72_example",   "prefix": "DCTdcParam",   "key_len": 3, "start_col": 3},
    # DCGEO columns (0-based): Id Name X Y Z TA RA1 RA2 L Res(9) ... Ofs(12)
    "residual":   {"dir": "DCGEO", "tpl": "DCGeomParam_e72_example",  "prefix": "DCGeomParam",  "key_len": 1, "start_col": 12},
    "resolution": {"dir": "DCGEO", "tpl": "DCGeomParam_e72_example",  "prefix": "DCGeomParam",  "key_len": 1, "start_col": 9},
}

RESI_CONV_MM = 0.010   # residual converged when every plane's |correction| is below this
RESI_IMPROVE = 0.8     # a residual pass improves when the largest correction shrinks below this fraction
RES_CONV = 0.05        # resolution converged when |p - 1| is below this
RES_STABLE_DP = 0.01   # a resolution pass does not change p when |dp| is below this

RES_GROUPS = {  # DCGEO Ids of each Res group
    "BLC1": range(1, 17), "BLC1a": range(1, 9), "BLC1b": range(9, 17),
    "BLC2": range(1001, 1017), "BLC2a": range(1001, 1009), "BLC2b": range(1009, 1017),
    "ALL": [*range(1, 17), *range(1001, 1017)],
}


def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("suffix", help="K or Pi")
    p.add_argument("ptype", choices=list(CATEGORIES))
    g = p.add_mutually_exclusive_group()
    g.add_argument("--bcout", action="store_true", help="BcOut (default)")
    g.add_argument("--bcin", action="store_true", help="BcIn")
    g.add_argument("--all", action="store_true", help="resolution only: one Res for BLC1+BLC2")
    p.add_argument("--separate", action="store_true", help="resolution only: one Res each for a and b")
    p.add_argument("--ftof", action="store_true", help="hdprm/hdphc: the FTOF detector set")
    a = p.parse_args()
    if a.suffix not in ("K", "Pi"):
        runner.fail(f"suffix should be K or Pi, not {a.suffix}")
    if a.ftof and a.ptype not in ("hdprm", "hdphc", "t0"):
        runner.fail("--ftof is for hdprm/hdphc only")
    if (a.all or a.separate) and a.ptype != "resolution":
        runner.fail("--all and --separate are for resolution only")
    a.chamber = "BcIn" if a.bcin else "BcOut"
    a.blc = "BLC1" if a.bcin else "BLC2"
    return a


# -- fit outputs -----------------------------------------------------------------
def fit_output(run, det, kind, suffix):
    name = "_".join(x for x in (f"run{run:05d}", det, kind, suffix) if x) + ".root"
    return config.OUTPUT_DIR / "root" / f"run{run:05d}" / name


def decode(run, tag):
    return config.DATA_DIR / f"run{run:05d}_{tag}.root"


def fresh_outputs(run, items, inputs, required=False):
    """{label: path} of items = [(label, det, kind, suffix)]. Stops if one is older than the newest input
    decode, or missing when required; a missing optional one is skipped with a warning."""
    t_input = max((p.resolve().stat().st_mtime for p in inputs if p.exists()), default=0.0)
    found, problems = {}, []
    for label, det, kind, suffix in items:
        path = fit_output(run, det, kind, suffix)
        if not path.exists():
            if required:
                problems.append(f"{label}: no fit output {path.name}")
            else:
                runner.warn(f"{label}: no fit output ({path.name}); its rows keep their current values")
        elif path.stat().st_mtime < t_input:
            problems.append(f"{label}: {path.name} is older than the decode (fit again)")
        else:
            found[label] = path
    if problems:
        runner.fail("nothing written:\n  " + "\n  ".join(problems))
    return found


def ok(label, n, note=""):
    print(f"  - {label}: [{colored('OK', 'green', attrs=['bold'])}] {n} entries found{note}")


# -- hodoscopes --------------------------------------------------------------------
def collect_hdprm(a, pf):
    dets = config.hodo_detectors("hdprm", a.ftof)
    found = fresh_outputs(a.run, [(d, d, "HDPRM", a.suffix) for d in dets], [decode(a.run, "Hodo")])
    data = {}
    for det, path in found.items():
        d = update_hdprm.make_dictdata(str(path), det)
        data.update(d)
        ok(f"{det:<6}", len(d))
    return data, []


def collect_t0(a, pf):
    found = fresh_outputs(a.run, [("T0_Offset", "T0_Offset", "", a.suffix)], [decode(a.run, "Hodo")])
    if "T0_Offset" not in found:
        return {}, []
    data = update_hdprm.make_t0_offset_dictdata(str(found["T0_Offset"]))
    ok("T0_Offset ", len(data))
    return data, []


def htof_inputs(run):
    return [decode(run, "Hodo"), decode(run, "TPCHelixHTOF")]


def collect_htofprm(a, pf):
    """HTOF TDC (all events, suffix-independent) and ADC (pion-tagged MIP) from HTOF_Calib."""
    found = fresh_outputs(a.run, [("HTOF TDC", "HTOF", "HDPRM", "All"), ("HTOF ADC", "HTOF", "ADC", a.suffix)],
                          htof_inputs(a.run))
    data = {}
    if "HTOF TDC" in found:
        d = update_hdprm.make_dictdata(str(found["HTOF TDC"]), "HTOF")
        data.update(d)
        ok("HTOF TDC", len(d))
    if "HTOF ADC" in found:
        path = str(found["HTOF ADC"])
        d = update_hdprm.make_htofprm_dictdata(path)
        data.update(d)
        ok("HTOF ADC", len(d))
        for seg, side, flag, mip in update_hdprm.htof_adc_quality_report(path):
            runner.warn(f"HTOF ADC seg{seg}{side}: fit-quality flag {flag}, MIP {mip:.1f} (written anyway)")
    return data, []


def phc_range(a, det):
    """Hand range of phc_conf.py; None = automatic choice in update_phc (HTOF: every segment)."""
    rng = phc_conf.limits_dict.get(f"{a.run:05d}_{a.suffix}_{det.lower()}")
    return [-np.inf, np.inf] if rng is None and det == "HTOF" else rng


def htof_phc_report(path):
    lim, off = update_hdprm.htof_phc_quality_report(path)
    for seg, side, p0, p1 in lim:
        runner.warn(f"HTOF seg{seg}{side}: PHC at fit limit (p0 {p0:.3f}, p1 {p1:.3f}), written anyway")
    for seg, n in off:
        print(colored(f"      HTOF seg{seg}: absolute TOF offset NOT applied (only {n} pion entries); "
                      "TOF of this segment is off", "red"))


def collect_phc(a, dets, inputs):
    found = fresh_outputs(a.run, [(d, d, "PHC", a.suffix) for d in dets], inputs)
    data = {}
    for det, path in found.items():
        rng = phc_range(a, det)
        try:
            d, not_used = update_phc.make_dictdata(str(path), det, good_ch_range=rng)
        except RuntimeError as e:
            runner.fail(f"{det}: {e}; nothing written")
        data.update(d)
        how = f"hand range {rng}" if rng is not None else "automatic"
        ok(f"{det:<6}", len(d), f" ({how}; mean of used segments for {not_used})")
        if det == "HTOF":
            htof_phc_report(str(path))
    return data, []


def collect_hdphc(a, pf):
    return collect_phc(a, config.hodo_detectors("hdphc", a.ftof), [decode(a.run, "Hodo")])


def collect_htofphc(a, pf):
    return collect_phc(a, ["HTOF"], htof_inputs(a.run))


# -- DC ----------------------------------------------------------------------------
def collect_dctdc(a, pf):
    path = fresh_outputs(a.run, [(a.blc, a.blc, "TDC", a.suffix)], [decode(a.run, a.chamber)], required=True)[a.blc]
    data = update_dctdc.make_dictdata(str(path))
    ok(f"{a.chamber} TDC", len(data))
    return data, []


def report_residual(corr, sigma, state):
    c = np.array(list(corr.values()))
    print(colored("\n>>> Convergence Report (Residual) <<<", "white", attrs=["bold"]))
    print(f"  - Correction (RMS):     {np.sqrt(np.mean(c ** 2)):.4f} mm")
    print(f"  - Maximum Correction:   {np.max(np.abs(c)):.4f} mm")
    print(f"  - Average Resolution:   {np.mean(list(sigma.values())):.4f} mm")
    status = {
        conv.CONVERGED: (f"[CONVERGED] (every |correction| < {RESI_CONV_MM * 1e3:.0f} um); not written", "green"),
        conv.SETTLED: (f"[SETTLED] max correction not shrinking over {conv.SETTLE_PASSES} passes; "
                       "not written, look at the residual PDF", "red"),
        conv.WRITTEN: ("[IN PROGRESS] written; decode again and run resi again", "yellow"),
    }[state]
    print(colored(f"  - Status: {status[0]}\n", status[1], attrs=["bold"]))


def collect_residual(a, pf):
    """New Ofs = current Ofs + correction (BLC_residual, relative to plane U1a)."""
    path = fresh_outputs(a.run, [(a.blc, a.blc, "residual", a.suffix)], [decode(a.run, a.chamber)],
                         required=True)[a.blc]
    fit = update_residual.make_dictdata(str(path))
    current = pf.get_data_dict(key_len=1)
    col = CATEGORIES["residual"]["start_col"] - 1
    missing = sorted(k for k in fit if k not in current)
    if not fit or missing:
        runner.fail(f"{path.name}: no planes, or Ids not in the DCGEO: {missing}")
    corr = {k: float(v[0]) for k, v in fit.items()}
    ok(f"{a.chamber} Residuals", len(corr), " (Additive Update)")
    max_corr = max(abs(c) for c in corr.values())
    state = conv.judge(a.run, a.blc, "resi", a.suffix, max_corr, max_corr < RESI_CONV_MM,
                       lambda prev, new: new > RESI_IMPROVE * prev)
    report_residual(corr, {k: float(v[1]) for k, v in fit.items()}, state)
    if state != conv.WRITTEN:
        return {}, [state]
    return {k: [float(current[k][col]) + c] for k, c in corr.items()}, [state]


def read_pull(path, need_ab):
    """pull_sigma, mean_chi2 (and pull_sigma_a/b) of a BLC_pull output."""
    import uproot
    branches = ["pull_sigma", "mean_chi2"] + (["pull_sigma_a", "pull_sigma_b"] if need_ab else [])
    f = uproot.open(str(path))
    tree = f["tree"] if "tree" in f else None
    missing = [b for b in branches if tree is None or b not in tree.keys()]
    if missing:
        runner.fail(f"{path.name}: no {missing} (rebuild BLC_pull and fit again)")
    arr = tree.arrays(branches, library="np")
    if len(arr["pull_sigma"]) < 1:
        runner.fail(f"{path.name}: empty tree")
    return {b: float(arr[b][0]) for b in branches}


def pull_groups(a):
    """[(group, p, mean_chi2)] of the Res groups chosen by --all / --separate."""
    blcs, tags = (["BLC1", "BLC2"], ["BcIn", "BcOut"]) if a.all else ([a.blc], [a.chamber])
    found = fresh_outputs(a.run, [(b, b, "pull", a.suffix) for b in blcs], [decode(a.run, t) for t in tags],
                          required=True)
    pulls = {b: read_pull(found[b], a.separate) for b in blcs}
    for b, v in pulls.items():
        ab = f" p_a={v['pull_sigma_a']:.4f} p_b={v['pull_sigma_b']:.4f}" if a.separate else ""
        runner.info(f"{b}: p={v['pull_sigma']:.4f}{ab} mean_chi2={v['mean_chi2']:.4f}")
    if a.separate:
        return [(f"{b}{s}", pulls[b][f"pull_sigma_{s}"], pulls[b]["mean_chi2"]) for b in blcs for s in "ab"]
    if a.all:
        return [("ALL", float(np.mean([v["pull_sigma"] for v in pulls.values()])),
                 float(np.mean([v["mean_chi2"] for v in pulls.values()])))]
    return [(a.blc, pulls[a.blc]["pull_sigma"], pulls[a.blc]["mean_chi2"])]


def report_resolution(group, p, chi2, state, r0, res_new):
    print(colored(f"\n>>> Convergence Report (Resolution / exclusive Pull) ({group}) <<<", "white", attrs=["bold"]))
    print(f"  - pull σ (p):           {p:.4f}  (target ≈ 1)")
    print(f"  - mean χ²_ν (monitor):  {chi2:.4f}")
    if r0 is not None:
        print(f"  - R0 / Res_new:         {r0:.6g} -> {res_new:.6g} mm")
    status = {
        conv.CONVERGED: (f"[CONVERGED] (|p-1|={abs(p - 1):.4f} < {RES_CONV}); not written", "green"),
        conv.SETTLED: (f"[SETTLED] p not changing over {conv.SETTLE_PASSES} passes but |p-1| >= {RES_CONV}; "
                       "not written, look at the pull PDF", "red"),
        conv.WRITTEN: (f"[IN PROGRESS] (|p-1|={abs(p - 1):.4f}) written; decode again and run reso again", "yellow"),
    }[state]
    print(colored(f"  - Status: {status[0]}\n", status[1], attrs=["bold"]))


def collect_resolution(a, pf):
    """Per group: judge p; if still moving, every layer gets Res_new = p * mean(current Res of the group)."""
    current = pf.get_data_dict(key_len=1)
    col = CATEGORIES["resolution"]["start_col"] - 1
    data, states = {}, []
    for group, p, chi2 in pull_groups(a):
        state = conv.judge(a.run, group, "reso", a.suffix, p, abs(p - 1.0) < RES_CONV,
                           lambda prev, new: abs(new - prev) < RES_STABLE_DP)
        states.append(state)
        r0 = res_new = None
        if state == conv.WRITTEN:
            keys = [str(i) for i in RES_GROUPS[group]]
            missing = [k for k in keys if k not in current]
            if missing:
                runner.fail(f"DCGEO has no Ids {missing} of {group}")
            r0 = float(np.mean([float(current[k][col]) for k in keys]))
            res_new = p * r0
            data.update({k: [f"{res_new:.6g}"] for k in keys})
        report_resolution(group, p, chi2, state, r0, res_new)
    return data, states


COLLECT = {
    "hdprm": collect_hdprm, "t0": collect_t0, "htofprm": collect_htofprm,
    "hdphc": collect_hdphc, "htofphc": collect_htofphc,
    "dctdc": collect_dctdc, "residual": collect_residual, "resolution": collect_resolution,
}


def new_param_source(cat, run, suffix):
    """Source for a run without its own param file: the existing param of the nearest run with the same
    suffix, else with the other suffix, else the template. Rows that are not written then keep that run's
    values instead of the template placeholders."""
    pdir = config.PARAM_DIR / cat["dir"] / config.SUB_DIR
    for s in (suffix, "Pi" if suffix == "K" else "K"):
        runs = {}
        for p in pdir.glob(f"{cat['prefix']}_run*_{s}"):
            m = re.fullmatch(rf"{cat['prefix']}_run(\d+)_{s}", p.name)
            if m:
                runs[int(m.group(1))] = p
        if runs:
            return runs[min(runs, key=lambda r: (abs(r - run), r))]
    return config.PARAM_DIR / cat["dir"] / cat["tpl"]


def main():
    a = parse_args()
    cat = CATEGORIES[a.ptype]
    target = config.PARAM_DIR / cat["dir"] / config.SUB_DIR / f"{cat['prefix']}_run{a.run:05d}_{a.suffix}"
    if not target.exists():
        source = new_param_source(cat, a.run, a.suffix)
        runner.warn(f"{target.name} does not exist; created from {source.name}")
        ParamFile.create_from_template(source, target)
    pf = ParamFile(target)
    print(colored(f"\n>>> Collecting results for {a.ptype} <<<", "cyan"))
    data, states = COLLECT[a.ptype](a, pf)
    if data:
        n = pf.update(data, key_len=cat["key_len"], start_col=cat["start_col"])
        pf.write()
        print(colored(f"\n[SUCCESS] Updated {n} entries total for Run {a.run} {a.suffix}", "green", attrs=["bold"]))
        print(f"File: {target.relative_to(config.ANALYZER_DIR)}\n")
    elif not states:
        runner.warn(f"no results for run {a.run} {a.suffix} {a.ptype}; nothing written")
    sys.exit(conv.exit_code(states))


if __name__ == "__main__":
    main()
