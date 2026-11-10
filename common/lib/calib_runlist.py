"""Helpers to call create_runlist and check Hodo TDC windows for calib chains."""
from pathlib import Path

from . import beam_run_condition as brc
from . import config, runner

SCRIPT_DIR = Path(__file__).resolve().parents[2] / "calibration" / "scripts"
CREATE_RUNLIST = SCRIPT_DIR / "create_runlist.py"

# Keys compared to the create_runlist USER template; if all match, window is unset.
HODO_TDC_KEYS = (
    "BHT_TDC", "BHT_TDC_Pi", "BHT_TDC_K",
    "BH2_TDC", "BAC_TDC", "HTOF_TDC", "KVC_TDC", "T1_TDC",
)


#_____________________________________________________________________________
def particle_prefix(kaon):
    return "K" if kaon else "Pi"


#_____________________________________________________________________________
def hodo_suffix(kaon, stage):
    """stage: 0 | hdprm | t0 | hdphc"""
    if stage == "0" or stage == 0:
        return "0"
    return f"{particle_prefix(kaon)}_{stage}"


#_____________________________________________________________________________
def dc_suffix(kaon, stage):
    """stage: tdc | drift | resi"""
    return f"{particle_prefix(kaon)}_{stage}"


#_____________________________________________________________________________
def d5_suffix(kaon):
    return particle_prefix(kaon)


#_____________________________________________________________________________
def resolve_d5_mom(run, explicit=None):
    """Momentum for create_runlist --d5. Prefer explicit, else beam_run_condition."""
    if explicit is not None:
        mom = float(explicit)
        if not (mom > 0):
            runner.fail("--mom must be > 0")
        return mom
    hit = brc.resolve(run)
    mom = hit.get("momentum_gev_c")
    if mom is None or not (float(mom) > 0):
        runner.fail(
            f"run{run:05d}: D5 create_runlist needs momentum "
            f"(table/UserParam or pass --mom)"
        )
    return float(mom)


#_____________________________________________________________________________
def user_param_path(run):
    for path in brc.user_param_paths(run):
        if path.is_file():
            return path
    return brc.user_param_paths(run)[0]


#_____________________________________________________________________________
def user_param_template_path():
    """Template that create_runlist copies for new USER files."""
    candidates = [
        config.PARAM_DIR / "USER" / "UserParam_e72_20251104",
        config.PARAM_DIR / "USER" / "example" / "UserParam_e72_example",
        config.PARAM_DIR / "USER" / "UserParam_e72_example",
    ]
    for path in candidates:
        if path.is_file():
            return path
    return candidates[0]


#_____________________________________________________________________________
def parse_tdc_windows(path):
    """Return {key: (min, max)} for lines like 'BHT_TDC  min  max'."""
    out = {}
    if path is None or not Path(path).is_file():
        return out
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        raw = line.split("#", 1)[0].strip()
        if not raw:
            continue
        parts = raw.split()
        if len(parts) < 3:
            continue
        key = parts[0]
        if not key.endswith("_TDC") and "_TDC_" not in key:
            # BHT_TDC_Pi / BHT_TDC_K included via _TDC_
            if not (key.startswith("BHT_TDC") or key.endswith("TDC")):
                continue
        try:
            lo = float(parts[1])
            hi = float(parts[2])
        except ValueError:
            continue
        out[key] = (lo, hi)
    return out


#_____________________________________________________________________________
def tdc_windows_ready(run):
    """Return (ready, user_path, detail).

    Ready when UserParam exists and at least one HODO_TDC_KEYS value differs
    from the create_runlist USER template (human has set windows for this run).
    """
    path = None
    for candidate in brc.user_param_paths(run):
        if candidate.is_file():
            path = candidate
            break
    if path is None:
        return False, user_param_path(run), "UserParam missing"
    tpl = user_param_template_path()
    if not tpl.is_file():
        # No template to compare: require keys present.
        got = parse_tdc_windows(path)
        missing = [k for k in HODO_TDC_KEYS if k not in got]
        if missing:
            return False, path, f"missing TDC keys {missing}"
        return True, path, "UserParam present (no template to compare)"
    user = parse_tdc_windows(path)
    ref = parse_tdc_windows(tpl)
    missing = [k for k in HODO_TDC_KEYS if k not in user]
    if missing:
        return False, path, f"missing TDC keys {missing}"
    same = 0
    differ = 0
    for key in HODO_TDC_KEYS:
        if key not in ref:
            continue
        if user.get(key) == ref.get(key):
            same += 1
        else:
            differ += 1
    if differ == 0:
        return False, path, f"TDC windows still match template {tpl.name}"
    return True, path, f"TDC windows differ from {tpl.name} ({differ} keys)"


#_____________________________________________________________________________
def announce_tdc_window(run):
    runner.warn(
        f"run{run:05d}: Hodo TDC window is not set yet (UserParam). "
        f"Calibration starts from UserParam / create_runlist."
    )
    runner.info(
        "Set windows, then re-run the hodo step:\n"
        f"  1) ~/tdcwin.sh {run}   # or calibration/scripts/set_tdc_window.py "
        f"<DATA_DIR>/run{run:05d}_Hodo.root BHT …\n"
        f"  2) ./ops.py decide {run} refit\n"
        f"  3) ./ops.py resume {run}   # or ./calib.py {run}"
    )


#_____________________________________________________________________________
def create_runlist_cmd(mode_flag, run, suffix, *, mom=None, dry_run=False):
    """Build shell command for create_runlist.py."""
    parts = ["python3", str(CREATE_RUNLIST), mode_flag, str(run), str(suffix)]
    if mom is not None:
        parts += ["--mom", str(mom)]
    if dry_run:
        parts.append("--dry-run")
    return " ".join(parts)


#_____________________________________________________________________________
def run_create_runlist(mode_flag, run, suffix, *, mom=None, dry_run=False):
    """Run create_runlist; fail on non-zero. Returns the command string."""
    if not CREATE_RUNLIST.is_file():
        runner.fail(f"create_runlist not found: {CREATE_RUNLIST}")
    cmd = create_runlist_cmd(mode_flag, run, suffix, mom=mom, dry_run=False)
    if dry_run:
        print(f"[dry-run] {cmd}")
        return cmd
    runner.step(f"create_runlist {mode_flag} {run} {suffix}")
    rc, _ = runner.sh(cmd, ok=(runner.EXIT_OK, runner.EXIT_FAIL))
    if rc != runner.EXIT_OK:
        runner.fail(f"create_runlist failed (rc={rc}): {cmd}")
    return cmd
