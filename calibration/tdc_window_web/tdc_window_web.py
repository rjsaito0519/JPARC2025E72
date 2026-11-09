#!/usr/bin/env python3
"""Set the Hodo TDC windows (USER param) by clicking in a browser; no X forwarding.

Usage: tdc_window_web.py [run] [--port 8765] [--trial]
  Starts a small web server on localhost; open the printed URL through an ssh tunnel
  (tdcwin_local.sh does both from the own PC). The run is chosen on the page (any run with
  DATA_DIR/runXXXXX_Hodo.root); the optional argument is only the run shown first.
  Click twice on a histogram for the window (BHT: three times, the third is the pi/K boundary),
  then press "write". The analyzer UserParam_runXXXXX is written (backup first, once per run and
  session); a missing one is made from the template, as set_tdc_window.py does.
  --trial writes copies in TRIAL_DIR instead.
"""
import argparse
import datetime
import gzip
import json
import re
import shutil
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np
import uproot

MYANALYSIS = Path(__file__).resolve().parents[2]
sys.path.append(str(MYANALYSIS))
from common.lib import config  # noqa: E402
from common.lib.param import ParamFile  # noqa: E402

HERE = Path(__file__).resolve().parent
BACKUP_BASE = Path("/group/had/sks/Users/sryuta/backup")
TRIAL_DIR = Path("/group/had/sks/Users/sryuta/tmp/tdc_window_web")
PLOTLY_JS = Path(np.__file__).parents[1] / "plotly" / "package_data" / "plotly.min.js"
USER_DIR = config.PARAM_DIR / "USER" / config.SUB_DIR
TEMPLATE = config.PARAM_DIR / "USER" / "UserParam_e72_example"
UD_DETS = ["BHT", "T0", "T1", "BH2", "HTOF", "CVC", "SAC3", "SFV"]
SUM_DETS = ["BAC", "KVC"]
TAB_ORDER = ["BHT", "T0", "BH2", "BAC", "HTOF", "KVC", "T1", "CVC", "SAC3", "SFV"]  # upstream first


# -- histograms -----------------------------------------------------------------
def hodo_root(run):
    return config.DATA_DIR / f"run{run:05d}_Hodo.root"


def available_runs():
    """[{run, decode}] for every DATA_DIR/runXXXXX_Hodo.root whose target exists."""
    out = []
    for p in sorted(config.DATA_DIR.glob("run[0-9][0-9][0-9][0-9][0-9]_Hodo.root")):
        if p.exists():
            out.append({"run": int(p.name[3:8]), "decode": p.resolve().name})
    return out


def sum_hists(f, keys):
    """(edges, summed counts) of the histograms in keys that exist in f; None if none exists."""
    edges, total = None, None
    for k in keys:
        if k not in f:
            continue
        counts, e = f[k].to_numpy()
        edges, total = (e, counts.copy()) if total is None else (edges, total + counts)
    return None if total is None else (edges, total)


def load_detector(f, det):
    """{"edges", "lines": {label: counts}} for one detector, summed as the old click tool does."""
    n = config.NUM_OF_CH[det]
    if det in UD_DETS:
        lines = {ud: sum_hists(f, [f"{det}_TDC_seg{i}{ud}" for i in range(n)]) for ud in ("U", "D")}
    elif det == "BAC":
        lines = {"sum": sum_hists(f, [f"BAC_TDC_seg{n}U"]) or sum_hists(f, [f"BAC_TDC_seg{i}U" for i in range(n)])}
    else:  # KVC: SUM channel of every segment
        lines = {"sum": sum_hists(f, [f"KVC_TDC_seg{i}S" for i in range(n)])}
    lines = {k: v for k, v in lines.items() if v is not None}
    if not lines:
        return None
    edges = next(iter(lines.values()))[0]
    return {"edges": edges.tolist(), "lines": {k: v[1].tolist() for k, v in lines.items()}}


def load_run(run):
    f = uproot.open(hodo_root(run))
    dets = {d: load_detector(f, d) for d in TAB_ORDER}
    return {d: v for d, v in dets.items() if v is not None}


# -- UserParam ------------------------------------------------------------------
def real_param(run):
    return USER_DIR / f"UserParam_run{run:05d}"


def target_param(run, real):
    return real_param(run) if real else TRIAL_DIR / real_param(run).name


def shown_param(run, real):
    """File whose windows are shown: the target, else the real one, else the template."""
    for p in (target_param(run, real), real_param(run), TEMPLATE):
        if p.exists():
            return p
    return None


def current_windows(path):
    if path is None:
        return {}
    data = ParamFile(path).get_data_dict()
    return {k: [float(v[0]), float(v[1])] for k, v in data.items() if k.endswith(("_TDC", "_TDC_Pi")) and len(v) >= 2}


def prepare_target(run, real, backups):
    """Make the target file if missing; back up the real file once per run and session."""
    target = target_param(run, real)
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(real_param(run) if real_param(run).exists() and not real else TEMPLATE, target)
    elif real and run not in backups:
        backup = BACKUP_BASE / f"userparam_tdcwin_{datetime.datetime.now():%Y%m%d_%H%M%S}"
        backup.mkdir(parents=True, exist_ok=True)
        shutil.copy2(target, backup / target.name)
        backups[run] = str(backup)
    return target


def window_values(det, lo, hi, sep):
    lo, hi = sorted([int(lo + 0.5), int(hi + 0.5)])
    values = {f"{det}_TDC": [lo, hi]}
    if sep is not None:
        s = int(sep + 0.5)
        values.update({f"{det}_TDC_Pi": [lo, s], f"{det}_TDC_K": [s, hi]})
    return values


def write_window(run, real, backups, req):
    target = prepare_target(run, real, backups)
    values = window_values(req["det"], req["lo"], req["hi"], req.get("sep"))
    pfile = ParamFile(target)
    n = pfile.update(values, key_len=1)
    pfile.write()
    return {"written": values, "target": str(target), "updated_rows": n, "expected_rows": len(values),
            "backup": backups.get(run)}


# -- server ---------------------------------------------------------------------
class App:
    def __init__(self, first_run, real):
        self.first_run, self.real = first_run, real
        self.cache, self.backups, self.lock = {}, {}, threading.Lock()

    def runs(self):
        return {"runs": available_runs(), "first": self.first_run, "real": self.real}

    def data(self, run):
        if not hodo_root(run).exists():
            return {"error": f"{hodo_root(run)} がありません"}
        with self.lock:
            if run not in self.cache:
                print(f"[load] run {run}: {hodo_root(run).resolve().name}")
                self.cache[run] = load_run(run)
        shown = shown_param(run, self.real)
        return {"run": run, "decode": str(hodo_root(run).resolve()), "target": str(target_param(run, self.real)),
                "shown": str(shown), "dets": self.cache[run], "current": current_windows(shown)}

    def write(self, req):
        with self.lock:
            res = write_window(int(req["run"]), self.real, self.backups, req)
        print(f"[write] run {req['run']}: {res['written']} -> {res['target']} "
              f"(rows {res['updated_rows']}/{res['expected_rows']})")
        return res


PLOTLY_RAW = PLOTLY_JS.read_bytes()
PLOTLY_GZ = gzip.compress(PLOTLY_RAW, 9)


def make_handler(app, page):
    class Handler(BaseHTTPRequestHandler):
        def send(self, body, ctype, cache=False):
            self.send_response(200)
            self.send_header("Content-Type", ctype)
            if "gzip" in self.headers.get("Accept-Encoding", ""):
                body = PLOTLY_GZ if body is PLOTLY_RAW else gzip.compress(body, 6)
                self.send_header("Content-Encoding", "gzip")
            if cache:
                self.send_header("Cache-Control", "max-age=2592000")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def send_json(self, obj):
            self.send(json.dumps(obj).encode(), "application/json")

        def do_GET(self):
            url = urlparse(self.path)
            if url.path == "/plotly.min.js":
                self.send(PLOTLY_RAW, "application/javascript", cache=True)
            elif url.path == "/runs":
                self.send_json(app.runs())
            elif url.path == "/data":
                run = parse_qs(url.query).get("run", [""])[0]
                self.send_json(app.data(int(run)) if run.isdigit() else {"error": f"run 番号が不正: {run}"})
            else:
                self.send(page.encode(), "text/html; charset=utf-8")

        def do_POST(self):
            self.send_json(app.write(json.loads(self.rfile.read(int(self.headers["Content-Length"])))))

        def log_message(self, *args):
            pass

    return Handler


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", nargs="?", default=None, help="run shown first (optional)")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--trial", action="store_true", help=f"write copies in {TRIAL_DIR} instead of the analyzer UserParam")
    args = ap.parse_args()
    first = int(re.sub(r"\D", "", args.run)) if args.run and re.search(r"\d", args.run) else None

    app = App(first, not args.trial)
    page = (HERE / "tdc_window_web.html").read_text()
    server = ThreadingHTTPServer(("127.0.0.1", args.port), make_handler(app, page))
    print(f"{len(available_runs())} runs with a Hodo root; writing to "
          f"{'trial copies in ' + str(TRIAL_DIR) if args.trial else 'the analyzer UserParam ' + str(USER_DIR)}")
    print(f"open http://localhost:{args.port}  (Ctrl-C to stop)")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("stopped")


if __name__ == "__main__":
    main()
