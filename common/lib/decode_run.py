"""Run myrun and verify decode success (stat DONE + fresh ROOT outputs).

myrun can exit 0 even when runs FAILED. Iterators must not treat process exit
alone as success. Mirrors the idea of iterate_tpc_calib's stat / mtime checks,
without polling (callers wait for myrun to finish, then verify).

Default is bsub (no --local). Pass local=True for work-server --local mode.
"""
import json
import time
from pathlib import Path

import yaml

from . import config, runner

MTIME_SLACK_SEC = 2.0
STAT_DIR = config.ANALYZER_DIR / "runmanager" / "stat"


#_____________________________________________________________________________
def decode_cmd(runlist, *, local=False):
    """Shell command: cd analyzer root and run myrun with an absolute runlist."""
    path = Path(runlist).expanduser().resolve()
    local_arg = " --local" if local else ""
    return f"cd {config.ANALYZER_DIR} && ./runmanager/myrun.py {path}{local_arg}"


#_____________________________________________________________________________
def stat_path(yml):
    return STAT_DIR / f"{Path(yml).stem}.json"


#_____________________________________________________________________________
def _workdir(doc):
    return Path(str(doc.get("WORKDIR", config.ANALYZER_DIR))).expanduser()


#_____________________________________________________________________________
def expected_runs_and_roots(yml):
    """Return (run_keys, absolute root paths) from the runlist RUN map."""
    path = Path(yml).expanduser().resolve()
    doc = yaml.safe_load(path.read_text())
    if not isinstance(doc, dict) or not isinstance(doc.get("RUN"), dict):
        runner.fail(f"runlist has no RUN map: {path}")
    wd = _workdir(doc)
    keys = []
    roots = []
    for key, entry in doc["RUN"].items():
        if not isinstance(entry, dict) or "root" not in entry:
            continue
        root = Path(str(entry["root"])).expanduser()
        if not root.is_absolute():
            root = wd / root
        keys.append(str(key))
        roots.append(root.resolve())
    if not keys:
        runner.fail(f"runlist has no RUN entries with root: {path}")
    return keys, roots


#_____________________________________________________________________________
def read_stat(path):
    if not path.is_file():
        return None
    raw = path.read_text(encoding="utf-8").strip()
    if not raw:
        return None
    try:
        info = json.loads(raw)
    except json.JSONDecodeError:
        return None
    return info if isinstance(info, dict) else None


#_____________________________________________________________________________
def analyze_stat(info, expected_keys):
    """Return (failed, unfinished, done_count)."""
    failed = []
    unfinished = []
    done_count = 0
    for key in expected_keys:
        item = info.get(key) if info else None
        if not isinstance(item, dict):
            unfinished.append(f"{key}:missing")
            continue
        st = str(item.get("stat", ""))
        if "FAIL" in st.upper():
            failed.append(f"{key}:{st}")
        elif st == "DONE":
            done_count += 1
        else:
            unfinished.append(f"{key}:{st or '?'}")
    return failed, unfinished, done_count


#_____________________________________________________________________________
def check_roots_fresh(roots, t0):
    """Return (missing, stale) absolute path strings."""
    cutoff = t0 - MTIME_SLACK_SEC
    missing = []
    stale = []
    for root in roots:
        if not root.is_file():
            missing.append(str(root))
            continue
        if root.stat().st_mtime < cutoff:
            stale.append(str(root))
    return missing, stale


#_____________________________________________________________________________
def verify_decode(yml, t0, myrun_rc=0):
    """Fail unless all expected runs are DONE and ROOT outputs are fresh."""
    keys, roots = expected_runs_and_roots(yml)
    sp = stat_path(yml)
    info = read_stat(sp)
    if info is None:
        runner.fail(f"decode stat missing or unreadable: {sp}")
    failed, unfinished, done_count = analyze_stat(info, keys)
    if failed or unfinished or done_count != len(keys):
        parts = []
        if failed:
            parts.append("FAILED " + ", ".join(failed))
        if unfinished:
            parts.append("unfinished " + ", ".join(unfinished))
        parts.append(f"done={done_count}/{len(keys)}")
        if myrun_rc != runner.EXIT_OK:
            parts.append(f"myrun_rc={myrun_rc}")
        runner.fail(f"decode not successful ({sp.name}): " + "; ".join(parts))
    missing, stale = check_roots_fresh(roots, t0)
    if missing or stale:
        parts = []
        if missing:
            parts.append("missing " + ", ".join(missing))
        if stale:
            parts.append("stale " + ", ".join(stale))
        runner.fail(f"decode outputs not fresh after myrun: " + "; ".join(parts))
    if myrun_rc != runner.EXIT_OK:
        # Stat and roots look fine but myrun reported failure — still stop.
        runner.fail(f"myrun exited {myrun_rc} despite DONE/fresh roots ({sp.name})")
    runner.info(f"decode OK: {sp.name} done={done_count}/{len(keys)}")


#_____________________________________________________________________________
def run_decode(runlist, *, dry_run=False, local=False):
    """Run myrun and verify stat DONE + fresh roots. dry_run: print only."""
    yml = Path(runlist).expanduser().resolve()
    cmd = decode_cmd(yml, local=local)
    if dry_run:
        print(f"[dry-run] {cmd}")
        return runner.EXIT_OK
    if not yml.is_file():
        runner.fail(f"runlist not found: {yml}")
    t0 = time.time()
    rc, _ = runner.sh(cmd, ok=(runner.EXIT_OK, runner.EXIT_FAIL))
    verify_decode(yml, t0, myrun_rc=rc)
    return runner.EXIT_OK
