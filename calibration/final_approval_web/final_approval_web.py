#!/usr/bin/env python3
"""Final calib QA approval Web UI (#25).

Usage: final_approval_web.py [run] [--port 8766]
  localhost only; open via SSH tunnel (final_approval_local.sh).
"""
from __future__ import annotations

import json
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

MYANALYSIS = Path(__file__).resolve().parents[2]
sys.path.append(str(MYANALYSIS))
from common.lib import (  # noqa: E402
    calib_approval,
    channel_flags,
    config,
    final_qa_index,
    qa_session,
)

HERE = Path(__file__).resolve().parent
HTML = HERE / "final_approval_web.html"
DISPATCH = MYANALYSIS / "calibration" / "scripts" / "calib_dispatch.py"
SESSION = lambda run: (
    Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / "calib_session.yml"
)


#_____________________________________________________________________________
def state_for(run: int) -> dict:
    path = calib_approval.dispatch_state_path(run)
    if not path.is_file():
        return {}
    try:
        st = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return st if isinstance(st, dict) else {}


#_____________________________________________________________________________
def available_runs() -> list[int]:
    root = Path(config.OUTPUT_DIR) / "calib"
    if not root.is_dir():
        return []
    runs = []
    for p in sorted(root.glob("run[0-9][0-9][0-9][0-9][0-9]")):
        if (p / "calib_dispatch_state.json").is_file():
            runs.append(int(p.name[3:8]))
    return runs


#_____________________________________________________________________________
def snapshot(run: int) -> dict:
    calib_approval.reconcile_stale_approval(run)
    st = state_for(run)
    subject = calib_approval.build_subject(run, st) if st else {}
    sha = calib_approval.subject_sha(subject) if subject else ""
    return {
        "run": run,
        "status": calib_approval.effective_status(st) if st else "none",
        "subject_sha": sha,
        "missing": list(subject.get("missing") or []) if subject else [],
        "hold": st.get("hold") if isinstance(st.get("hold"), dict) else None,
        "unverified": list(st.get("unverified") or []),
        "open_flags": channel_flags.open_needs_review(run),
        "pdf_index": (
            final_qa_index.index_from_subject(run, subject) if subject
            else {"run": run, "n_pdf": 0, "pdfs": [], "source": "subject_sealed"}
        ),
        "has_session": SESSION(run).is_file(),
    }


#_____________________________________________________________________________
def run_approve_or_revoke(
    run: int,
    op: str,
    reason: str = "",
    *,
    subject_sha: str = "",
) -> tuple[bool, str]:
    yml = SESSION(run)
    if not yml.is_file():
        return False, f"no session YAML: {yml}"
    import subprocess
    cmd = [sys.executable, str(DISPATCH), op, str(yml), "--via", "web"]
    if op == "approve" and subject_sha:
        cmd += ["--subject-sha", str(subject_sha)]
    if op == "revoke" and reason:
        cmd += ["--reason", reason]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"{type(exc).__name__}: {exc}"
    text = ((result.stdout or "") + (result.stderr or "")).strip()
    if len(text) > 1500:
        text = "…\n" + text[-1500:]
    return result.returncode == 0, text or f"exit={result.returncode}"


#_____________________________________________________________________________
class Handler(BaseHTTPRequestHandler):
    first_run = None

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
        if u.path == "/runs":
            self._json(200, {
                "runs": available_runs(),
                "first": self.first_run,
            })
            return
        if u.path == "/state":
            qs = parse_qs(u.query)
            try:
                run = int((qs.get("run") or [0])[0])
            except ValueError:
                self._json(400, {"error": "bad run"})
                return
            if run < 1:
                self._json(400, {"error": "bad run"})
                return
            self._json(200, snapshot(run))
            return
        if u.path == "/pdf":
            qs = parse_qs(u.query)
            name = (qs.get("name") or [""])[0]
            exec_id = (qs.get("exec_id") or [""])[0]
            try:
                run = int((qs.get("run") or [0])[0])
            except ValueError:
                self._send(400, b"bad run", "text/plain")
                return
            if not name or not exec_id or run < 1:
                self._send(400, b"bad pdf query", "text/plain")
                return
            man_root = qa_session.manifest_dir(run).resolve()
            man_path = man_root / f"{exec_id}.json"
            sealed = (man_root / exec_id / name).resolve()
            try:
                sealed.relative_to(man_root)
            except ValueError:
                self._send(403, b"forbidden", "text/plain")
                return
            if not man_path.is_file() or not sealed.is_file():
                self._send(404, b"not found", "text/plain")
                return
            if sealed.suffix.lower() != ".pdf":
                self._send(404, b"not found", "text/plain")
                return
            try:
                manifest = qa_session.load_manifest(man_path)
            except (OSError, ValueError, TypeError):
                self._send(404, b"not found", "text/plain")
                return
            ok, _msg, paths = qa_session.resolve_sealed_outputs(
                manifest, run=run,
            )
            if not ok or sealed not in {p.resolve() for p in paths}:
                self._send(409, b"sealed sha mismatch", "text/plain")
                return
            self._send(200, sealed.read_bytes(), "application/pdf")
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
        try:
            run = int(body.get("run") or 0)
        except (TypeError, ValueError):
            self._json(400, {"error": "bad run"})
            return
        if run < 1:
            self._json(400, {"error": "bad run"})
            return
        if u.path == "/approve":
            subject_sha = str(body.get("subject_sha") or "").strip()
            if not subject_sha:
                self._json(400, {
                    "ok": False,
                    "error": "subject_sha required",
                    "state": snapshot(run),
                })
                return
            ok, msg = run_approve_or_revoke(
                run, "approve", subject_sha=subject_sha,
            )
            self._json(200 if ok else 400, {"ok": ok, "message": msg, "state": snapshot(run)})
            return
        if u.path == "/revoke":
            ok, msg = run_approve_or_revoke(
                run, "revoke", str(body.get("reason") or "web revoke"),
            )
            self._json(200 if ok else 400, {"ok": ok, "message": msg, "state": snapshot(run)})
            return
        self._json(404, {"error": "not found"})


#_____________________________________________________________________________
def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", nargs="?", type=int, default=None)
    parser.add_argument("--port", type=int, default=8766)
    a = parser.parse_args()
    Handler.first_run = a.run
    host = "127.0.0.1"
    httpd = ThreadingHTTPServer((host, a.port), Handler)
    print(f"final_approval_web on {host}:{a.port}; open http://{host}:{a.port}/")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nbye")


if __name__ == "__main__":
    main()
