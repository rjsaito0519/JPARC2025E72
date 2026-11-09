"""Shared pieces of the calibration front ends (run_hodo.py, run_dcprm.py, run_d5.py).

Every front end takes `<run...> <mode> [--kaon] [--debug] [--force]` and processes the runs one by
one. Exit codes for a driver:
  EXIT_OK            calibration succeeded (non-iterative write, or iterative converged)
  EXIT_FAIL          stopped; nothing more was written after the failure
  EXIT_NOT_CONVERGED params written, another decode + iteration is needed
  EXIT_CHECK         not converged / settled, or QA check failed; nothing useful written, a person looks
  EXIT_SKIP          finished without applying calibration (no fit results, or --debug)
"""
import argparse
import subprocess
import sys
from pathlib import Path

from termcolor import colored

from . import config
from .input_check import check_input
from .param import backup_file  # noqa: F401  (re-exported for the front ends)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = PROJECT_ROOT / "bin"
SCRIPT_DIR = PROJECT_ROOT / "calibration" / "scripts"

EXIT_OK, EXIT_FAIL, EXIT_NOT_CONVERGED, EXIT_CHECK, EXIT_SKIP = 0, 1, 2, 3, 4

# keep our lines in order with the output of the commands we run (logs, pipes, batch jobs)
sys.stdout.reconfigure(line_buffering=True)


def info(msg):
    print(colored(f"[INFO] {msg}", "green"))


def warn(msg):
    print(colored(f"[Warning] {msg}", "yellow"))


def step(msg):
    print(colored(f">>> {msg}", "cyan"))


def fail(msg):
    print(colored(f"[Error] {msg}", "red"))
    sys.exit(EXIT_FAIL)


def sh(cmd, capture=False, check=True, ok=(EXIT_OK,)):
    """Run a shell command; (return code, output if capture). check=True stops unless the code is in ok."""
    print(colored(f"Running: {cmd}", "cyan"))
    res = subprocess.run(cmd, shell=True, capture_output=capture, text=True)
    out = (res.stdout or "") + (res.stderr or "")
    if check and res.returncode not in ok:
        if out:
            print(out)
        fail(f"command failed with exit code {res.returncode}")
    return res.returncode, out


class Parser(argparse.ArgumentParser):
    """ArgumentParser whose usage errors exit with EXIT_FAIL (argparse's own 2 means EXIT_NOT_CONVERGED here)."""

    def __init__(self, description):
        super().__init__(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)

    def error(self, message):
        self.print_usage(sys.stderr)
        fail(message)


def base_parser(description, modes):
    """Parser with the arguments every front end shares."""
    p = Parser(description)
    p.add_argument("runs", type=int, nargs="+", help="run number(s), processed one by one")
    p.add_argument("mode", choices=modes)
    p.add_argument("--kaon", action="store_true", help="K (param suffix _K) instead of pi (_Pi)")
    p.add_argument("--debug", action="store_true", help="fit and make PDFs only; no param is written")
    p.add_argument("--force", action="store_true",
                   help="run even if the input decode is older than its params or is of the other particle")
    return p


def suffix_of(args):
    return "K" if args.kaon else "Pi"


def require_bins(*names):
    """Paths of the fit binaries; stops before anything runs if one is not built."""
    paths = [BIN_DIR / n for n in names]
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        fail(f"not built: {', '.join(missing)} (run build.sh)")
    return paths


def checked_root(path, suffix, force):
    """path if it exists and passes check_input; stops otherwise."""
    if not Path(path).exists():
        fail(f"{path} not found")
    if not check_input(path, suffix, force):
        sys.exit(EXIT_FAIL)
    info(f"input: {path}")
    return Path(path)


def input_root(run, tag, suffix, force):
    """DATA_DIR/runXXXXX_<tag>.root (the link create_runlist.py makes), checked."""
    path = config.DATA_DIR / f"run{run:05d}_{tag}.root"
    if not path.exists():
        fail(f"{path} not found (create_runlist.py for {tag}, then decode)")
    return checked_root(path, suffix, force)


def update_param(run, suffix, ptype, *opts, debug=False):
    """update_param.py <run> <suffix> <ptype> [opts]; skipped with --debug (EXIT_SKIP).

    Child exit codes accepted without stopping: 0, 2, 3, 4.
    """
    if debug:
        warn(f"--debug: {ptype} params not written")
        return EXIT_SKIP
    step(f"update params ({ptype})")
    rc, _ = sh(" ".join(["python3", str(SCRIPT_DIR / "update_param.py"), str(run), suffix, ptype, *opts]),
               ok=(EXIT_OK, EXIT_NOT_CONVERGED, EXIT_CHECK, EXIT_SKIP))
    return rc


def worst(codes):
    """Exit code for several groups or runs: CHECK > NOT_CONVERGED > SKIP > OK."""
    order = (EXIT_OK, EXIT_SKIP, EXIT_NOT_CONVERGED, EXIT_CHECK)
    return max(codes, key=order.index, default=EXIT_OK)
