#!/usr/bin/env python3
"""Interactive Hodo HDPRM fit Web (#23). BH2 first; other dets / modes share API.

Usage: hodo_fit_web.py [run] [--port 8767] [--suffix Pi]
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import uproot

MYANALYSIS = Path(__file__).resolve().parents[2]
sys.path.append(str(MYANALYSIS))
from common.lib import channel_flags, config, fit_overrides, runner  # noqa: E402

HERE = Path(__file__).resolve().parent
HTML = HERE / "hodo_fit_web.html"
DEFAULT_DET = "BH2"
DEFAULT_MODE = "hdprm"


#_____________________________________________________________________________
def hodo_root(run: int, suffix: str) -> Path:
    return Path(config.DATA_DIR) / f"run{int(run):05d}_Hodo.root"


#_____________________________________________________________________________
def fit_root(run: int, det: str, suffix: str) -> Path:
    return (
        Path(config.OUTPUT_DIR) / "root" / f"run{int(run):05d}"
        / f"run{int(run):05d}_{det}_HDPRM_{suffix}.root"
    )


#_____________________________________________________________________________
def hist_name(det: str, seg: int, side: str, kind: str, suffix: str) -> str:
    return f"{det}_{kind}_seg{int(seg)}{side.upper()}_{suffix}"


#_____________________________________________________________________________
def load_hist(run: int, det: str, seg: int, side: str, kind: str, suffix: str):
    path = hodo_root(run, suffix)
    if not path.is_file():
        return None, f"missing {path.name}"
    name = hist_name(det, seg, side, kind, suffix)
    with uproot.open(path) as f:
        if name not in f:
            return None, f"missing hist {name}"
        counts, edges = f[name].to_numpy()
    return {
        "name": name,
        "edges": edges.tolist(),
        "counts": counts.tolist(),
    }, ""


#_____________________________________________________________________________
def load_fit_vals(run: int, det: str, seg: int, suffix: str):
    path = fit_root(run, det, suffix)
    if not path.is_file():
        return None
    with uproot.open(path) as f:
        tree = f["tree"].arrays(library="np")
    for i, ch in enumerate(tree["ch"]):
        if int(ch) != int(seg):
            continue
        return {
            "adc_ped": [float(tree["adc_p0_val"][i][0]), float(tree["adc_p0_val"][i][1])],
            "adc_mip": [float(tree["adc_p1_val"][i][0]), float(tree["adc_p1_val"][i][1])],
            "tdc": [float(tree["tdc_p0_val"][i][0]), float(tree["tdc_p0_val"][i][1])],
        }
    return None


#_____________________________________________________________________________
def n_seg(det: str) -> int:
    return int(config.NUM_OF_CH.get(det.upper()) or 8)


#_____________________________________________________________________________
def param_key(det: str, seg: int, aort: int, ud: int) -> str:
    cid = config.DETECTOR_I_LIST[det]
    return f"{cid}-0-{int(seg)}-{int(aort)}-{int(ud)}"


#_____________________________________________________________________________
def refit_bh2(run: int, suffix: str) -> tuple[bool, str]:
    flat = fit_overrides.write_flat_mip_left(run)
    os.environ["E72_FIT_OVERRIDE"] = str(flat)
    try:
        (exe,) = runner.require_bins("BH2_HDPRM")
    except SystemExit as exc:
        return False, str(exc)
    root = runner.input_root(run, "Hodo", suffix, force=False)
    env = os.environ.copy()
    env["E72_FIT_OVERRIDE"] = str(flat)
    try:
        result = subprocess.run(
            [str(exe), str(root), suffix],
            capture_output=True, text=True, env=env, timeout=600,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"{type(exc).__name__}: {exc}"
    text = ((result.stdout or "") + (result.stderr or "")).strip()
    if len(text) > 2000:
        text = "…\n" + text[-2000:]
    return result.returncode == 0, text or f"exit={result.returncode}"


#_____________________________________________________________________________
def adopt_hdprm(run: int, suffix: str) -> tuple[bool, str]:
    rc = runner.update_param(run, suffix, "hdprm")
    return rc == 0, f"update_param hdprm rc={rc}"


#_____________________________________________________________________________
class Handler(BaseHTTPRequestHandler):
    first_run = None
    suffix = "Pi"

    def log_message(self, fmt, *args):
        sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    def _send(self, code, body, content_type="application/json"):
        data = body if isinstance(body, (bytes, bytearray)) else str(body).encode()
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _json(self, code, obj):
        self._send(code, json.dumps(obj, ensure_ascii=False), "application/json")

    def do_GET(self):
        u = urlparse(self.path)
        if u.path in ("/", "/index.html"):
            self._send(200, HTML.read_bytes(), "text/html; charset=utf-8")
            return
        if u.path == "/meta":
            self._json(200, {
                "first": self.first_run,
                "suffix": self.suffix,
                "dets": ["BH2", "CVC", "T1", "SAC3"],
                "modes": ["hdprm", "t0", "hdphc", "htofprm"],
                "default_det": DEFAULT_DET,
                "default_mode": DEFAULT_MODE,
            })
            return
        if u.path == "/channel":
            qs = parse_qs(u.query)
            try:
                run = int((qs.get("run") or [0])[0])
                seg = int((qs.get("seg") or [0])[0])
            except ValueError:
                self._json(400, {"error": "bad run/seg"})
                return
            det = (qs.get("det") or [DEFAULT_DET])[0]
            side = (qs.get("side") or ["U"])[0]
            kind = (qs.get("kind") or ["ADC"])[0]
            mode = (qs.get("mode") or [DEFAULT_MODE])[0]
            suffix = (qs.get("suffix") or [self.suffix])[0]
            hist, err = load_hist(run, det, seg, side, kind, suffix)
            doc = fit_overrides.load(run)
            key = fit_overrides.channel_key(det, seg, side)
            ov = (
                ((doc.get("fit") or {}).get(det.upper()) or {}).get(mode) or {}
            ).get(key) or {}
            locks = fit_overrides.calib_lock_keys(run, det, mode=mode)
            self._json(200, {
                "run": run,
                "det": det,
                "seg": seg,
                "side": side,
                "mode": mode,
                "n_seg": n_seg(det),
                "hist": hist,
                "hist_error": err,
                "fit": load_fit_vals(run, det, seg, suffix),
                "override": ov,
                "fit_key": key,
                "param_keys": {
                    "adc": param_key(det, seg, 0, 0 if side.upper() == "U" else 1),
                    "tdc": param_key(det, seg, 1, 0 if side.upper() == "U" else 1),
                },
                "calib_locks": sorted(locks),
            })
            return
        self._send(404, b"not found", "text/plain")

    def do_POST(self):
        u = urlparse(self.path)
        length = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(length) if length else b"{}"
        try:
            body = json.loads(raw.decode() or "{}")
        except json.JSONDecodeError:
            self._json(400, {"error": "bad json"})
            return
        run = int(body.get("run") or 0)
        if run < 1:
            self._json(400, {"error": "bad run"})
            return
        det = str(body.get("det") or DEFAULT_DET)
        mode = str(body.get("mode") or DEFAULT_MODE)
        suffix = str(body.get("suffix") or self.suffix)

        if u.path == "/override":
            seg = int(body.get("seg") or 0)
            side = str(body.get("side") or "U")
            mip_left = float(body["mip_left"])
            doc = fit_overrides.set_mip_left(run, det, seg, side, mip_left, mode=mode)
            self._json(200, {"ok": True, "doc": doc})
            return
        if u.path == "/refit":
            if det.upper() != "BH2" and mode == "hdprm":
                # Other HDPRM dets also work via env; BH2 is the primary smoke path.
                pass
            if mode != "hdprm":
                self._json(400, {"ok": False, "message": f"mode {mode} refit not wired yet"})
                return
            flat = fit_overrides.write_flat_mip_left(run)
            os.environ["E72_FIT_OVERRIDE"] = str(flat)
            try:
                (exe,) = runner.require_bins(f"{det}_HDPRM")
            except SystemExit as exc:
                self._json(400, {"ok": False, "message": str(exc)})
                return
            root = runner.input_root(run, "Hodo", suffix, force=False)
            env = os.environ.copy()
            env["E72_FIT_OVERRIDE"] = str(flat)
            try:
                result = subprocess.run(
                    [str(exe), str(root), suffix],
                    capture_output=True, text=True, env=env, timeout=600,
                )
            except (OSError, subprocess.TimeoutExpired) as exc:
                self._json(400, {"ok": False, "message": f"{type(exc).__name__}"})
                return
            text = ((result.stdout or "") + (result.stderr or "")).strip()
            self._json(200, {
                "ok": result.returncode == 0,
                "message": text[-2000:] if text else f"exit={result.returncode}",
            })
            return
        if u.path == "/adopt":
            ok, msg = adopt_hdprm(run, suffix)
            self._json(200 if ok else 400, {"ok": ok, "message": msg})
            return
        if u.path == "/lock":
            key = str(body.get("param_key") or "")
            locked = bool(body.get("locked"))
            if not key:
                self._json(400, {"error": "param_key required"})
                return
            doc = fit_overrides.set_calib_lock(run, det, key, locked, mode=mode)
            self._json(200, {"ok": True, "doc": doc})
            return
        if u.path == "/flag":
            channel = str(body.get("channel") or "")
            reason = str(body.get("reason") or "web flag")
            channel_flags.append_flag(
                run,
                family="hodo",
                mode=mode,
                det=det,
                channel=channel,
                source="human",
                reason=reason,
                status="needs_review",
            )
            self._json(200, {"ok": True})
            return
        self._json(404, {"error": "not found"})


#_____________________________________________________________________________
def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", nargs="?", type=int, default=None)
    parser.add_argument("--port", type=int, default=8767)
    parser.add_argument("--suffix", default="Pi", choices=("Pi", "K"))
    a = parser.parse_args()
    Handler.first_run = a.run
    Handler.suffix = a.suffix
    host = "127.0.0.1"
    httpd = ThreadingHTTPServer((host, a.port), Handler)
    print(f"hodo_fit_web on {host}:{a.port}; open http://{host}:{a.port}/")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nbye")


if __name__ == "__main__":
    main()
