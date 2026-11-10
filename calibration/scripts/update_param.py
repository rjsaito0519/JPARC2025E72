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
exit codes (lib/runner.py): 0 calib OK (converged or non-iterative write), 1 fail, 2 written (iterate),
3 settled / check, 4 skip (no results).
"""
import math
import re
import sys
from pathlib import Path

import numpy as np
from termcolor import colored

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import calib_anomaly, config, runner  # noqa: E402
from common.lib.param import ParamFile  # noqa: E402
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
def _hdprm_skip_anomaly(a, det, path):
    """Refuse write when BAC/KVC skip ratio is high; return reason or None."""
    import uproot
    tree = uproot.open(str(path))["tree"].arrays(library="np")
    if det == "BAC":
        bad, total, labels = calib_anomaly.bac_skip_counts(tree)
        family = "hdprm/BAC"
    elif det == "KVC":
        bad, total, labels = calib_anomaly.kvc_skip_counts(tree)
        family = "hdprm/KVC"
    else:
        return None
    reason = calib_anomaly.decide_ratio(
        family, bad, total, calib_anomaly.HDPRM_SKIP_RATIO, labels,
    )
    if reason:
        calib_anomaly.remember_check_reason(a.run, reason)
        calib_anomaly.dump_candidates(
            a.run,
            calib_anomaly.candidates_payload("hodo", "hdprm", det, labels, reason),
        )
    return reason


def collect_hdprm(a, pf):
    dets = config.hodo_detectors("hdprm", a.ftof)
    found = fresh_outputs(a.run, [(d, d, "HDPRM", a.suffix) for d in dets], [decode(a.run, "Hodo")])
    data = {}
    for det, path in found.items():
        if _hdprm_skip_anomaly(a, det, path):
            return {}, [conv.SETTLED]
        d = update_hdprm.make_dictdata(str(path), det, run=a.run)
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
        flags = update_hdprm.htof_adc_quality_report(path)
        for seg, side, flag, mip in flags:
            runner.warn(f"HTOF ADC seg{seg}{side}: fit-quality flag {flag}, MIP {mip:.1f} (written anyway)")
        n_flag, labels = calib_anomaly.htof_adc_flag_counts(flags)
        # 34 segs x 3 sides is the usual HTOF layout; fall back to flagged-only total.
        total_sides = max(34 * 3, n_flag)
        reason = calib_anomaly.decide_ratio(
            "htofprm/ADC", n_flag, total_sides, calib_anomaly.HTOF_FLAG_RATIO, labels,
        )
        if reason:
            calib_anomaly.remember_check_reason(a.run, reason)
            calib_anomaly.dump_candidates(
                a.run,
                calib_anomaly.candidates_payload(
                    "hodo", "htofprm", "HTOF", labels, reason,
                ),
            )
            return {}, [conv.SETTLED]
        d = update_hdprm.make_htofprm_dictdata(path)
        d = update_hdprm._drop_locked(a.run, "HTOF", d, mode="htofprm")
        data.update(d)
        ok("HTOF ADC", len(d))
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
    from common.lib import phc_range_gate

    found = fresh_outputs(a.run, [(d, d, "PHC", a.suffix) for d in dets], inputs)
    data = {}
    proposals = []
    for det, path in found.items():
        hand_rng = phc_conf.limits_dict.get(
            f"{a.run:05d}_{a.suffix}_{det.lower()}"
        )
        rng = phc_range(a, det)
        try:
            d, not_used = update_phc.make_dictdata(str(path), det, good_ch_range=rng)
        except RuntimeError as e:
            runner.fail(f"{det}: {e}; nothing written")
        # Hand range: outside segs are intentional; do not treat as PHC anomaly.
        if hand_rng is not None:
            data.update(d)
            ok(
                f"{det:<6}", len(d),
                f" (hand range {hand_rng}; mean of used segments for {not_used})",
            )
            if det == "HTOF":
                htof_phc_report(str(path))
            continue
        import uproot
        tree = uproot.open(str(path))["tree"].arrays(library="np")
        n_seg = len(tree["ch"])
        bad, total, labels = calib_anomaly.phc_bad_counts(not_used, n_seg)
        reason = calib_anomaly.decide_ratio(
            f"phc/{det}", bad, total, calib_anomaly.PHC_BAD_RATIO, labels,
        )
        if reason:
            ch_list = [int(c) for c in tree["ch"]]
            proposed = phc_range_gate.propose_from_good(ch_list, not_used)
            if proposed is None:
                calib_anomaly.remember_check_reason(a.run, reason)
                calib_anomaly.dump_candidates(
                    a.run,
                    calib_anomaly.candidates_payload(
                        "hodo", "hdphc", det, labels, reason,
                    ),
                )
                return {}, [conv.SETTLED]
            first, last = proposed
            proposals.append({
                "det": det,
                "first": first,
                "last": last,
                "bad": bad,
                "total": total,
                "labels": labels[:32],
            })
            continue
        data.update(d)
        ok(f"{det:<6}", len(d), f" (automatic; mean of used segments for {not_used})")
        if det == "HTOF":
            htof_phc_report(str(path))
    if proposals:
        payload = {
            "run": int(a.run),
            "suffix": str(a.suffix),
            "proposals": proposals,
        }
        phc_range_gate.write_proposal(a.run, payload)
        ok_slack, slack_msg = phc_range_gate.notify_phc_range_slack(
            a.run, proposals,
        )
        if not ok_slack:
            runner.warn(f"PHC range Slack note failed: {slack_msg}")
        dets_s = ",".join(str(p["det"]) for p in proposals)
        ranges_s = "; ".join(
            f"{p['det']}[{p['first']},{p['last']}]" for p in proposals
        )
        calib_anomaly.remember_check_reason(
            a.run,
            f"phc_range_proposal: {dets_s} -> {ranges_s}",
        )
        for p in proposals:
            calib_anomaly.dump_candidates(
                a.run,
                calib_anomaly.candidates_payload(
                    "hodo", "hdphc", p["det"], p.get("labels") or [],
                    f"phc_range_proposal:{p['det']}[{p['first']},{p['last']}]",
                ),
            )
        # Nothing written until Accept / Custom fills phc_conf and refits.
        return {}, [conv.SETTLED]
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
    """New Ofs = current Ofs + correction (BLC_residual, relative to plane U1a).

    Returns (data, states, pending_records). Journal pending/written is recorded around the write.
    """
    decode_path = decode(a.run, a.chamber)
    path = fresh_outputs(a.run, [(a.blc, a.blc, "residual", a.suffix)], [decode_path],
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
    hist = conv.history_path(a.run, a.blc, "resi", a.suffix)
    apply_key = conv.make_apply_key(conv.input_token(decode_path), conv.corr_token(corr))
    state = conv.decide_at(hist, a.blc, max_corr, max_corr < RESI_CONV_MM,
                           lambda prev, new: new > RESI_IMPROVE * prev)
    report_residual(corr, {k: float(v[1]) for k, v in fit.items()}, state)
    if state == conv.WRITTEN:
        st = conv.apply_status(hist, a.blc, apply_key, pf.path)
        if st == conv.APPLY_DONE:
            runner.warn(f"{a.blc} resi: same decode+correction already applied; not writing again "
                        "(decode and refit before another residual pass)")
            return {}, [conv.SKIPPED], []
        if st == conv.APPLY_MISMATCH:
            runner.fail(
                f"{a.blc} resi: param content matches neither before nor after hash for this apply "
                f"(apply={apply_key}); inspect {hist} and {pf.path}"
            )
    if state in (conv.CONVERGED, conv.SETTLED):
        conv.record(hist, a.blc, max_corr, state, conv.details_with_apply(apply_key))
        if state == conv.SETTLED:
            calib_anomaly.remember_check_reason(
                a.run,
                f"residual SETTLED {a.blc}: max_corr={max_corr:.4f} mm "
                f"(not shrinking over {conv.SETTLE_PASSES} passes)",
            )
        return {}, [state], []
    if state != conv.WRITTEN:
        return {}, [state], []
    pending = [{"path": hist, "group": a.blc, "value": max_corr, "state": conv.WRITTEN,
                "apply_key": apply_key}]
    return {k: [float(current[k][col]) + c] for k, c in corr.items()}, [state], pending


def read_pull(path, need_ab):
    """pull_sigma, mean_chi2 (and pull_sigma_a/b) of a BLC_pull output.

    Optional pull_fit_ok / pull_fit_ok_a / pull_fit_ok_b are included when present
    (older ROOT files without them remain valid).
    """
    import uproot
    required = ["pull_sigma", "mean_chi2"] + (["pull_sigma_a", "pull_sigma_b"] if need_ab else [])
    optional = ["pull_fit_ok"] + (["pull_fit_ok_a", "pull_fit_ok_b"] if need_ab else [])
    f = uproot.open(str(path))
    tree = f["tree"] if "tree" in f else None
    missing = [b for b in required if tree is None or b not in tree.keys()]
    if missing:
        runner.fail(f"{path.name}: no {missing} (rebuild BLC_pull and fit again)")
    have_opt = [b for b in optional if b in tree.keys()]
    arr = tree.arrays(required + have_opt, library="np")
    if len(arr["pull_sigma"]) < 1:
        runner.fail(f"{path.name}: empty tree")
    out = {b: float(arr[b][0]) for b in required}
    for b in have_opt:
        out[b] = int(arr[b][0])
    return out


def _pull_fit_ok_combined(pulls, key="pull_fit_ok"):
    """None if no file has the flag; else 1 only when every present flag is non-zero."""
    flags = [v[key] for v in pulls.values() if key in v]
    if not flags:
        return None
    return 1 if all(f != 0 for f in flags) else 0


def pull_groups(a):
    """[(group, p, mean_chi2, fit_ok)] of the Res groups chosen by --all / --separate.

    fit_ok is None when the ROOT has no pull_fit_ok branch (legacy); else 0/1.
    """
    blcs, tags = (["BLC1", "BLC2"], ["BcIn", "BcOut"]) if a.all else ([a.blc], [a.chamber])
    found = fresh_outputs(a.run, [(b, b, "pull", a.suffix) for b in blcs], [decode(a.run, t) for t in tags],
                          required=True)
    pulls = {b: read_pull(found[b], a.separate) for b in blcs}
    for b, v in pulls.items():
        ab = f" p_a={v['pull_sigma_a']:.4f} p_b={v['pull_sigma_b']:.4f}" if a.separate else ""
        runner.info(f"{b}: p={v['pull_sigma']:.4f}{ab} mean_chi2={v['mean_chi2']:.4f}")
    if a.separate:
        return [(f"{b}{s}", pulls[b][f"pull_sigma_{s}"], pulls[b]["mean_chi2"],
                 pulls[b].get(f"pull_fit_ok_{s}"))
                for b in blcs for s in "ab"]
    if a.all:
        return [("ALL", float(np.mean([v["pull_sigma"] for v in pulls.values()])),
                 float(np.mean([v["mean_chi2"] for v in pulls.values()])),
                 _pull_fit_ok_combined(pulls))]
    return [(a.blc, pulls[a.blc]["pull_sigma"], pulls[a.blc]["mean_chi2"],
             pulls[a.blc].get("pull_fit_ok"))]


def pull_width_ok(p, fit_ok):
    """True when p can be used to scale DCGEO Res."""
    if fit_ok == 0:
        return False
    return math.isfinite(p) and p > 0


def report_resolution(group, p, chi2, state, r0, res_new, invalid=False):
    print(colored(f"\n>>> Convergence Report (Resolution / exclusive Pull) ({group}) <<<", "white", attrs=["bold"]))
    p_txt = f"{p:.4f}" if math.isfinite(p) else str(p)
    print(f"  - pull σ (p):           {p_txt}  (target ≈ 1)")
    chi2_txt = f"{chi2:.4f}" if math.isfinite(chi2) else str(chi2)
    print(f"  - mean χ²_ν (monitor):  {chi2_txt}")
    if r0 is not None:
        print(f"  - R0 / Res_new:         {r0:.6g} -> {res_new:.6g} mm")
    if invalid:
        print(colored("  - Status: [CHECK] invalid pull σ "
                      "(not finite, <=0, or pull_fit_ok=0); DCGEO not written\n",
                      "red", attrs=["bold"]))
        return
    status = {
        conv.CONVERGED: (f"[CONVERGED] (|p-1|={abs(p - 1):.4f} < {RES_CONV}); not written", "green"),
        conv.SETTLED: (f"[SETTLED] p not changing over {conv.SETTLE_PASSES} passes but |p-1| >= {RES_CONV}; "
                       "not written, look at the pull PDF", "red"),
        conv.WRITTEN: (f"[IN PROGRESS] (|p-1|={abs(p - 1):.4f}) written; decode again and run reso again", "yellow"),
    }[state]
    print(colored(f"  - Status: {status[0]}\n", status[1], attrs=["bold"]))


def resolution_apply_key(hist, group, decode_tok, p):
    """Identify a measured scale, independently of the Res values it updates.

    Old journals hashed the resulting Res, so they cannot safely identify a
    retry after a write. Require a new decode for those inputs during migration.
    """
    prefix = conv.make_apply_key(decode_tok, "")
    if hist.exists():
        for line in hist.read_text().splitlines():
            fields = line.split("\t")
            if len(fields) < 5 or fields[1] != group or fields[3] not in (conv.PENDING, conv.WRITTEN):
                continue
            for detail in fields[4].split():
                if detail.startswith("apply=" + prefix) and not detail.startswith("apply=" + prefix + "reso-v2:"):
                    runner.fail(f"{group} reso: legacy apply journal for this decode; "
                                "decode again before applying resolution (do not delete the journal)")
    corr = {"group": group, "ids": list(RES_GROUPS[group]), "p": p}
    return conv.make_apply_key(decode_tok, "reso-v2:" + conv.corr_token(corr))


def collect_resolution(a, pf):
    """Per group: decide p; if still moving, every layer gets Res_new = p * mean(current Res of the group).

    Returns (data, states, pending_records). Journal pending/written is recorded around the write.
    """
    current = pf.get_data_dict(key_len=1)
    col = CATEGORIES["resolution"]["start_col"] - 1
    data, states, pending = {}, [], []
    tags = ["BcIn", "BcOut"] if a.all else [a.chamber]
    decode_tok = "+".join(conv.input_token(decode(a.run, t)) for t in tags)
    for group, p, chi2, fit_ok in pull_groups(a):
        if not pull_width_ok(p, fit_ok):
            states.append(conv.SETTLED)
            report_resolution(group, p, chi2, conv.SETTLED, None, None, invalid=True)
            calib_anomaly.remember_check_reason(
                a.run,
                f"resolution CHECK {group}: invalid pull σ "
                f"(p={p}, fit_ok={fit_ok})",
            )
            continue
        hist = conv.history_path(a.run, group, "reso", a.suffix)
        state = conv.decide_at(hist, group, p, abs(p - 1.0) < RES_CONV,
                               lambda prev, new: abs(new - prev) < RES_STABLE_DP)
        r0 = res_new = None
        keys = []
        if state == conv.WRITTEN:
            keys = [str(i) for i in RES_GROUPS[group]]
            missing = [k for k in keys if k not in current]
            if missing:
                runner.fail(f"DCGEO has no Ids {missing} of {group}")
            r0 = float(np.mean([float(current[k][col]) for k in keys]))
            res_new = p * r0
        apply_key = resolution_apply_key(hist, group, decode_tok, p)
        report_resolution(group, p, chi2, state, r0, res_new)
        if state == conv.WRITTEN:
            st = conv.apply_status(hist, group, apply_key, pf.path)
            if st == conv.APPLY_DONE:
                runner.warn(f"{group} reso: same decode+correction already applied; not writing again "
                            "(decode and refit before another resolution pass)")
                states.append(conv.SKIPPED)
                continue
            if st == conv.APPLY_MISMATCH:
                runner.fail(
                    f"{group} reso: param content matches neither before nor after hash for this apply "
                    f"(apply={apply_key}); inspect {hist} and {pf.path}"
                )
        if state == conv.WRITTEN:
            data.update({k: [f"{res_new:.6g}"] for k in keys})
            pending.append({"path": hist, "group": group, "value": p, "state": conv.WRITTEN,
                            "apply_key": apply_key})
        elif state in (conv.CONVERGED, conv.SETTLED):
            conv.record(hist, group, p, state, conv.details_with_apply(apply_key))
            if state == conv.SETTLED:
                calib_anomaly.remember_check_reason(
                    a.run,
                    f"resolution SETTLED {group}: p={p:.4f} "
                    f"(not changing over {conv.SETTLE_PASSES} passes)",
                )
        states.append(state)
    return data, states, pending


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
    collected = COLLECT[a.ptype](a, pf)
    data, states = collected[0], collected[1]
    pending = collected[2] if len(collected) > 2 else []
    if data:
        missing, duplicates = pf.key_issues(data, key_len=cat["key_len"])
        if missing or duplicates:
            parts = []
            if missing:
                parts.append(f"missing keys: {missing}")
            if duplicates:
                parts.append(f"duplicate keys: {duplicates}")
            runner.fail("nothing written (" + "; ".join(parts) + ")")
        before = conv.file_sha256(target)
        n = pf.update(data, key_len=cat["key_len"], start_col=cat["start_col"])
        if n != len(data):
            runner.fail(f"nothing written (updated {n} rows, expected {len(data)})")
        after = conv.sha256_hex(pf.serialized_text())
        for rec in pending:
            key = rec.get("apply_key") or rec.get("fit_token")
            conv.record(
                rec["path"], rec["group"], rec["value"], conv.PENDING,
                conv.details_with_apply(key, before=before, after=after),
            )
        pf.write()
        after_disk = conv.file_sha256(target)
        for rec in pending:
            key = rec.get("apply_key") or rec.get("fit_token")
            conv.record(
                rec["path"], rec["group"], rec["value"], conv.WRITTEN,
                conv.details_with_apply(key, before=before, after=after_disk),
            )
        print(colored(f"\n[SUCCESS] Updated {n} entries total for Run {a.run} {a.suffix}", "green", attrs=["bold"]))
        print(f"File: {target.relative_to(config.ANALYZER_DIR)}\n")
    elif not states:
        runner.warn(f"no results for run {a.run} {a.suffix} {a.ptype}; nothing written")
        sys.exit(runner.EXIT_SKIP)
    sys.exit(conv.exit_code(states))


if __name__ == "__main__":
    main()
