"""Checks on a decoded ROOT before a calibration front end (run_dcprm.py, run_hodo.py) updates params.

- stale input (warn + y/N): a calibration param used by the decode (HDPRM, HDPHC, DCGEO, DCTDC,
  DCDRFT in the conf that produced the ROOT) is newer than the ROOT, i.e. the params changed and the
  run was not decoded again. Updates relative to the current values (DCGEO Ofs / Res) would be
  applied twice.
- particle mismatch (stop): the input link points to a decode of the other particle (_Pi for K,
  _K for Pi).

The conf is found from the runmanager runlists (the RUN entry whose root is the input file).
"""
import sys
from pathlib import Path

import yaml
from termcolor import colored

from . import config

CALIB_PARAM_KEYS = ("HDPRM", "HDPHC", "DCGEO", "DCTDC", "DCDRFT")
RUNLIST_DIR = config.ANALYZER_DIR / "runmanager" / "runlist"
_OTHER_SUFFIX = {"K": "Pi", "Pi": "K"}


def _workdir(doc):
    return Path(str(doc.get("WORKDIR", config.ANALYZER_DIR))).expanduser()


def find_conf_for_root(root_path):
    """Conf (absolute Path) of the runlist entry that writes root_path, or None."""
    target = Path(root_path).resolve()
    for yml in sorted(RUNLIST_DIR.glob("**/*.yml")):
        try:
            doc = yaml.safe_load(yml.read_text())
        except Exception:
            continue
        if not isinstance(doc, dict) or not isinstance(doc.get("RUN"), dict):
            continue
        wd = _workdir(doc)
        for entry in doc["RUN"].values():
            if not isinstance(entry, dict) or "root" not in entry or "conf" not in entry:
                continue
            root = Path(str(entry["root"])).expanduser()
            if not root.is_absolute():
                root = wd / root
            if root.resolve() == target:
                conf = Path(str(entry["conf"])).expanduser()
                return conf if conf.is_absolute() else wd / conf
    return None


def read_conf_params(conf_path, keys=CALIB_PARAM_KEYS):
    """{key: absolute Path} for the given param keys of an analyzer conf."""
    params = {}
    for line in Path(conf_path).read_text().splitlines():
        f = line.split()
        if len(f) >= 2 and f[0].rstrip(":") in keys:
            p = Path(f[1])
            params[f[0].rstrip(":")] = p if p.is_absolute() else config.ANALYZER_DIR / p
    return params


def _ask(question):
    """y/N prompt; False when there is no terminal to ask."""
    if not sys.stdin.isatty():
        print(colored("        (no terminal to ask: not continuing; use --force to run anyway)", "yellow"))
        return False
    try:
        return input(colored(f"        {question} [y/N]: ", "yellow")).strip().lower() in ("y", "yes")
    except EOFError:
        return False


def check_input(root_path, suffix, force=False, keys=CALIB_PARAM_KEYS):
    """Return True if root_path is fine to calibrate with.

    A decode of the other particle stops (unless force). Params newer than the decode only warn and
    ask y/N (the DCGEO of BLC1 and BLC2 is one file, so a BLC1 update also dates the BcOut decode).
    With force=True everything is reported as a warning and True is returned.
    """
    real = Path(root_path).resolve()

    other = _OTHER_SUFFIX.get(suffix)
    if other and other in real.stem.split("_"):
        msg = f"input is a {other} decode ({real.name}) but suffix {suffix} was requested"
        if not force:
            print(colored(f"[Error] {msg}\n        (use --force to run anyway)", "red"))
            return False
        print(colored(f"[Warning] {msg}", "yellow"))

    conf = find_conf_for_root(real)
    if conf is None or not conf.exists():
        print(colored(f"[Warning] no runlist entry writes {real}; staleness of the params not checked", "yellow"))
        return True
    t_root = real.stat().st_mtime
    newer = [f"{k}: {p}" for k, p in read_conf_params(conf, keys).items()
             if p.exists() and p.stat().st_mtime > t_root + 1.0]
    if not newer:
        return True
    print(colored(f"[Warning] params changed after the decode {real.name} (conf {conf.name}):\n      "
                  + "\n      ".join(newer), "yellow"))
    print(colored("        Updates relative to the current values (DCGEO Ofs, Res) would be applied twice "
                  "if this decode was already used.", "yellow"))
    if force:
        return True
    return _ask("Continue with this decode?")
