#!/usr/bin/env python3
"""BH2 MIP tuning over the same localhost SSH tunnel as the TDC window editor."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import re
import secrets
import shutil
import subprocess
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
import uproot

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from common.lib import config, fit_overrides  # noqa: E402
from common.lib.calib_state_io import locked_path, write_json  # noqa: E402

HERE = Path(__file__).resolve().parent
BINARY = PROJECT / "bin" / "BH2_MIP_PREVIEW"
PLOTLY = Path(np.__file__).parents[1] / "plotly/package_data/plotly.min.js"


def read_json(path, default=None):
    if not path.exists():
        return copy.deepcopy(default)
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


class App:
    def __init__(self, run, particle="Pi"):
        self.run, self.particle = int(run), particle
        self.base = Path(config.OUTPUT_DIR) / "calib" / f"run{run:05d}"
        self.work = Path(config.WORK_DIR) / "mip_web" / f"run{run:05d}" / particle
        self.work.mkdir(parents=True, exist_ok=True)
        self.staging_path = self.base / f"mip_web_{particle}.json"
        self.root = Path(config.DATA_DIR) / f"run{run:05d}_Hodo.root"
        self.lock = threading.RLock()
        self.preview_lock = threading.Lock()
        self.cache = {}
        self.token = secrets.token_urlsafe(32)

    def fingerprint(self):
        resolved = self.root.resolve(strict=True)
        stat = resolved.stat()
        param = Path(config.PARAM_DIR)
        sources = [fit_overrides.overrides_path(self.run),
                   param / "HDPRM" / config.SUB_DIR / f"HodoParam_run{self.run:05d}_{self.particle}",
                   param / "HDPHC" / config.SUB_DIR / f"HodoPHCParam_run{self.run:05d}_{self.particle}",
                   param / "USER" / config.SUB_DIR / f"UserParam_run{self.run:05d}", BINARY,
                   PROJECT / "common/include/params.h"]
        fit_root = Path(config.OUTPUT_DIR) / "root" / f"run{self.run:05d}" / f"run{self.run:05d}_BH2_HDPRM_{self.particle}.root"
        fit_stat = fit_root.stat() if fit_root.is_file() else None
        record = [str(resolved), stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns,
                  [(str(path), digest(path)) for path in sources],
                  [fit_stat.st_size, fit_stat.st_mtime_ns] if fit_stat else None]
        return hashlib.sha256(json.dumps(record).encode()).hexdigest()

    def channel(self, key):
        match = re.fullmatch(r"bh2-(\d+)-(u|d)", str(key))
        if not match or not 0 <= int(match[1]) < config.NUM_OF_CH["BH2"]:
            raise ValueError("BH2 channel is invalid")
        return int(match[1]), match[2].upper()

    def document(self):
        current = self.fingerprint()
        doc = read_json(self.staging_path, {"fingerprint": current, "staged": {}, "receipts": {}, "revision": 0})
        job = doc.get("job") or {}
        if job.get("status") == "running":
            state = read_json(Path(job["state_path"]), {})
            hold = state.get("hold") or {}
            if state.get("status") == "held" and hold.get("id") != job["hold_id"]:
                job["status"] = "done" if hold.get("step") == "hodo_hdprm" and hold.get("rc") == 0 else "failed"
                job["message"] = "再実行完了・Slackの新しいQAを確認してください" if job["status"] == "done" else "再実行が確認待ちで停止しました。Slackとログを確認してください"
            elif state.get("status") in ("failed", "stopped", "awaiting_final_approval", "approved") and hold.get("id") != job["hold_id"]:
                job["status"] = "done" if state.get("status") == "awaiting_final_approval" and state.get("last_rc", state.get("rc")) == 0 else "failed"
                job["message"] = "校正段階が終了しました。SlackでQAを確認してください" if job["status"] == "done" else "再実行が停止しました。Slackとログを確認してください"
            if (job["status"] == "running" and state.get("status") in ("held", "stopped")
                    and hold.get("id") == job["hold_id"]
                    and time.time() - job.get("started_at", time.time()) > 120):
                job.update(status="failed", message="設定は保存済みですが、再実行の開始を確認できません。Slackと再開ログを確認してください")
            if job["status"] != "running":
                doc["fingerprint"] = current
                doc["staged"], doc["receipts"] = {}, {}
                doc["revision"] = doc.get("revision", 0) + 1
                self.cache.clear()
                write_json(self.staging_path, doc)
        if doc["fingerprint"] != current and job.get("status") != "running":
            had_staged = bool(doc.get("staged"))
            doc.update(fingerprint=current, staged={}, receipts={}, revision=doc.get("revision", 0) + 1)
            doc["notice"] = "元のデータまたは設定が変わったため、以前の仮保存を破棄しました" if had_staged else "元のデータまたは設定が変わったため、結果を更新しました"
            self.cache.clear()
            write_json(self.staging_path, doc)
        return doc

    def gate(self):
        for path in [self.base / "lanes/hodo/calib_dispatch_state.json", self.base / "calib_dispatch_state.json"]:
            state = read_json(path, {})
            hold = state.get("hold") or {}
            if state.get("status") == "held" and hold.get("step") == "hodo_hdprm" and not hold.get("decision"):
                import yaml
                session_path = self.base / "calib_session.yml"
                session = yaml.safe_load(session_path.read_text()) or {}
                if session.get("auto"):
                    raise ValueError("自動進行のセッションではWriteできません。承認待ちを有効にしてください")
                from common.lib import calib_approval
                _, suffix = calib_approval.particle_suffix(state)
                if suffix != self.particle:
                    raise ValueError("表示中の粒子と校正セッションが異なります")
                return path, state, hold
        raise ValueError("WriteはHodo ADC/TDC（hodo_hdprm）の確認待ちから行ってください。現在は閲覧・調整・仮保存のみ可能です")

    def public(self, doc):
        try:
            _, _, hold = self.gate()
            gate = {"hold_id": hold["id"], "reason": ""}
        except (ValueError, OSError) as error:
            gate = {"hold_id": None, "reason": str(error)}
        return {"run": self.run, "particle": self.particle, "fingerprint": doc["fingerprint"],
                "revision": doc["revision"], "staged": doc["staged"], "notice": doc.get("notice", ""),
                "job": doc.get("job"), "gate": gate}

    def state(self):
        with self.lock, locked_path(self.staging_path):
            return self.public(self.document())

    def require_current(self, request, doc):
        if (doc.get("job") or {}).get("status") == "running":
            raise ValueError("再実行中は閲覧・編集・Writeを停止しています")
        if request.get("fingerprint") != doc["fingerprint"]:
            raise ValueError("元のデータや設定が変わりました。一覧を読み直してください")

    def execute_preview(self, key, settings=None, curve_only=False):
        segment, side = self.channel(key)
        identifier = secrets.token_hex(12)
        output = self.work / f"{identifier}.json"
        settings_path = self.work / f"{identifier}.txt"
        command = [str(BINARY), str(self.root.resolve()), self.particle, str(segment), side, str(output)]
        if settings is not None:
            settings_path.write_text(fit_overrides.mip_fit_line(key, self.particle, settings) + "\n")
            # The preview uses the same particle-qualified key as the production reader.
            if curve_only:
                command = [str(BINARY), "--curve", str(settings_path), f"{key}:{self.particle}", str(output)]
            else:
                command += [str(settings_path)]
        environment = dict(os.environ)
        environment["E72_FIT_OVERRIDE"] = str(fit_overrides.flat_override_path(self.run))
        try:
            with self.preview_lock:
                result = subprocess.run(command, env=environment, capture_output=True, text=True, timeout=90)
            if result.returncode:
                raise RuntimeError("Fit処理に失敗: " + (result.stderr or result.stdout)[-1000:])
            return read_json(output)
        finally:
            output.unlink(missing_ok=True)
            settings_path.unlink(missing_ok=True)

    def overview(self, request):
        with self.lock, locked_path(self.staging_path):
            doc = self.document()
            self.require_current(request, doc)
            offset = int(request.get("page", 0)) * 12
            keys = [f"bh2-{segment}-{side}" for segment in range(config.NUM_OF_CH["BH2"]) for side in ("u", "d")]
            if offset < 0 or offset >= len(keys):
                raise ValueError("page is invalid")
            with uproot.open(self.root) as source:
                channels = []
                for key in keys[offset:offset+12]:
                    segment, side = self.channel(key)
                    name = f"BH2_ADC_seg{segment}{side}_{self.particle}"
                    counts, edges = source[name].to_numpy()
                    if not np.all(np.isfinite(counts)):
                        raise ValueError(f"nonfinite histogram: {name}")
                    counts = counts.reshape(-1, 2).sum(axis=1)
                    edges = edges[::2]
                    channels.append({"key": key, "label": f"BH2 seg{segment} {side}",
                                     "x": ((edges[:-1] + edges[1:])/2).tolist(), "counts": counts.tolist(),
                                     "automatic": self.automatic(key, doc["fingerprint"])})
            return {**self.public(doc), "channels": channels, "pages": (len(keys)+11)//12}

    def automatic(self, key, fingerprint):
        cache_key = (fingerprint, key)
        if cache_key not in self.cache:
            fit_root = Path(config.OUTPUT_DIR) / "root" / f"run{self.run:05d}" / f"run{self.run:05d}_BH2_HDPRM_{self.particle}.root"
            snapshot = None
            if fit_root.is_file() and fit_root.stat().st_mtime_ns >= self.root.stat().st_mtime_ns:
                with uproot.open(fit_root) as output:
                    if "mip_fit" in output:
                        segment, side = self.channel(key)
                        columns = output["mip_fit"].arrays(library="np")
                        for index, (row_seg, row_side) in enumerate(zip(columns["seg"], columns["side"])):
                            if row_seg == segment and row_side == (0 if side == "U" else 1):
                                parameters = columns["parameters"][index].tolist()
                                bounds = columns["range"][index].tolist()
                                if bounds[0] < bounds[1] and parameters[0] > 0 and parameters[2] > 0:
                                    function = str(columns["function"][index])
                                    snapshot = {"parameters": parameters, "range": bounds, "function": function,
                                                "success": bool(columns["success"][index]), "status": int(columns["status"][index]),
                                                "covariance": int(columns["covariance"][index]), "ndf": int(columns["ndf"][index]),
                                                "origin": "前回実行のFit結果" if int(columns["status"][index]) >= 0 else "低統計の自動推定値（Fitなし）を初期値として表示"}
                                    curve_settings = {"parameters": parameters, "range": bounds,
                                                      "function": function, "range_mode": "fixed"}
                                    snapshot.update(self.execute_preview(key, curve_settings, curve_only=True))
                                break
            self.cache[cache_key] = snapshot or {**self.execute_preview(key), "origin": "自動Fitを現在のデータから再計算（校正値は更新しません）"}
        if self.cache[cache_key].get("status") == -1:
            self.cache[cache_key]["origin"] = "低統計の自動推定値（Fitなし）を初期値として表示"
        return copy.deepcopy(self.cache[cache_key])

    def edit(self, request):
        with self.lock, locked_path(self.staging_path):
            doc = self.document()
            self.require_current(request, doc)
            key = request["key"]
            self.channel(key)
            result = self.automatic(key, doc["fingerprint"])
            entry = fit_overrides.load(self.run).get("fit", {}).get("BH2", {}).get("hdprm", {}).get(key, {})
            committed = entry.get("mip_fit", {}).get(self.particle)
            automatic = {"parameters": result["parameters"], "range": result["range"],
                         "function": "auto", "range_mode": "auto"}
            initial = doc["staged"].get(key, {}).get("settings") or committed or automatic
            return {**self.public(doc), "automatic": result, "initial": initial,
                    "baseline": committed or automatic, "automatic_note": result["origin"]}

    def fit(self, request, curve_only=False):
        with self.lock, locked_path(self.staging_path):
            doc = self.document()
            self.require_current(request, doc)
            key = request["key"]
            settings = fit_overrides.validate_mip_fit(request["settings"])
            if not 0 <= settings["range"][0] < settings["range"][1] <= 4096:
                raise ValueError("Fit範囲はADCヒストグラム内（0–4096）で指定してください")
            if curve_only:
                effective = dict(settings)
                if effective["function"] == "auto":
                    effective["function"] = request.get("auto_function", "gaussian")
                return self.execute_preview(key, effective, curve_only=True)
            result = self.execute_preview(key, settings)
            self.require_current(request, self.document())
            if result["success"]:
                settings["parameters"] = result["parameters"]
                settings["range"] = result["range"]
                receipt = secrets.token_hex(16)
                doc["receipts"][receipt] = {"key": key, "settings": settings, "result": result}
                # Keep only the most recent successful fit for this channel.
                doc["receipts"] = {token: row for token, row in doc["receipts"].items()
                                   if row["key"] != key or token == receipt}
                write_json(self.staging_path, doc)
                result["receipt"] = receipt
                result["settings"] = settings
            return result

    def stage(self, request):
        with self.lock, locked_path(self.staging_path):
            doc = self.document()
            self.require_current(request, doc)
            key = request["key"]
            self.channel(key)
            if request.get("remove"):
                doc["staged"].pop(key, None)
            else:
                receipt = doc["receipts"].get(request.get("receipt"))
                if not receipt or receipt["key"] != key:
                    raise ValueError("成功したFitを実行してから仮保存してください")
                doc["staged"][key] = copy.deepcopy(receipt)
            doc["revision"] += 1
            write_json(self.staging_path, doc)
            return self.public(doc)

    def review(self, request):
        with self.lock, locked_path(self.staging_path):
            doc = self.document()
            self.require_current(request, doc)
            rows = []
            for key, row in doc["staged"].items():
                entry = fit_overrides.load(self.run).get("fit", {}).get("BH2", {}).get("hdprm", {}).get(key, {})
                before = entry.get("mip_fit", {}).get(self.particle)
                if before is None:
                    automatic = self.automatic(key, doc["fingerprint"])
                    before = {"parameters": automatic["parameters"], "range": automatic["range"],
                              "function": "auto", "range_mode": "auto"}
                after = row["settings"]
                changes = []
                for field in ("function", "range_mode", "range", "parameters"):
                    if before[field] != after[field]:
                        if isinstance(after[field], list):
                            for index, (old, new) in enumerate(zip(before[field], after[field])):
                                if old != new:
                                    label = f"p{index}" if field == "parameters" else ("下限", "上限")[index]
                                    changes.append([label, old, new])
                        else:
                            changes.append([field, before[field], after[field]])
                rows.append({"key": key, "changes": changes})
            return {**self.public(doc), "rows": rows}

    def write(self, request):
        # Keep settings and the replay receipt under the same cross-session lock.
        with self.lock, locked_path(self.staging_path), locked_path(fit_overrides.overrides_path(self.run)):
            doc = self.document()
            self.require_current(request, doc)
            if request.get("revision") != doc["revision"] or not doc["staged"]:
                raise ValueError("仮保存が変更されました。確認画面を開き直してください")
            state_path, state, hold = self.gate()
            with locked_path(state_path):
                latest = read_json(state_path, {})
                if latest != state:
                    raise ValueError("承認待ちの状態が変わりました。確認画面を開き直してください")
                if request.get("hold_id") != hold["id"]:
                    raise ValueError("承認待ちの対象が変わりました。確認画面を開き直してください")
                overrides = read_json(fit_overrides.overrides_path(self.run), fit_overrides.empty_doc())
                channels = overrides.setdefault("fit", {}).setdefault("BH2", {}).setdefault("hdprm", {})
                for key, row in doc["staged"].items():
                    channels.setdefault(key, {}).setdefault("mip_fit", {})[self.particle] = fit_overrides.validate_mip_fit(row["settings"])
                target = fit_overrides.overrides_path(self.run)
                if target.exists():
                    backup = self.work / (time.strftime("overrides_%Y%m%dT%H%M%S_") + secrets.token_hex(4) + ".json")
                    shutil.copy2(target, backup)
                fit_overrides.save(self.run, overrides)
                doc["job"] = {"status": "running", "hold_id": hold["id"], "state_path": str(state_path),
                              "started_at": time.time(), "message": "HDPRMを再実行中。完了後にSlackの新しいQAで確認待ちになります"}
                write_json(self.staging_path, doc)
                # Refit uses the existing lane's decode mode and always returns to its gate.
                from common.lib import calib_ops
                try:
                    accepted, message = calib_ops.run_dispatch("retry", self.run, hold_id=hold["id"], background_resume=True)
                    if not accepted:
                        raise RuntimeError(message)
                except BaseException as error:
                    doc["job"].update(status="failed", message=f"設定は保存済み。再実行の開始に失敗: {error}")
                    doc["fingerprint"] = self.fingerprint()
                    write_json(self.staging_path, doc)
                return self.public(doc)


def handler(app):
    class Handler(BaseHTTPRequestHandler):
        def send(self, payload, content_type="application/json", status=200):
            body = payload if isinstance(payload, bytes) else json.dumps(payload, ensure_ascii=False, allow_nan=False).encode()
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def valid_host(self):
            host = self.headers.get("Host", "").split(":")[0]
            if host not in ("localhost", "127.0.0.1"):
                self.send({"error": "localhost access only"}, status=403)
                return False
            return True

        def do_GET(self):
            if not self.valid_host():
                return
            try:
                if urlparse(self.path).path == "/plotly.min.js":
                    self.send(PLOTLY.read_bytes(), "application/javascript")
                elif urlparse(self.path).path == "/state":
                    self.send(app.state())
                elif urlparse(self.path).path == "/":
                    page = (HERE / "mip_web.html").read_text().replace("__TOKEN__", app.token)
                    self.send(page.encode(), "text/html; charset=utf-8")
                else:
                    self.send({"error": "not found"}, status=404)
            except Exception as error:
                self.send({"error": str(error)}, status=409)

        def do_POST(self):
            if not self.valid_host():
                return
            try:
                if self.headers.get("X-MIP-Token") != app.token:
                    self.send({"error": "invalid session"}, status=403)
                    return
                size = int(self.headers.get("Content-Length", "0"))
                if size < 1 or size > 65536:
                    raise ValueError("invalid request size")
                request = json.loads(self.rfile.read(size))
                methods = {"/overview": app.overview, "/edit": app.edit, "/fit": app.fit,
                           "/curve": lambda req: app.fit(req, curve_only=True),
                           "/stage": app.stage, "/review": app.review, "/write": app.write}
                method = methods.get(urlparse(self.path).path)
                if method is None:
                    self.send({"error": "not found"}, status=404)
                    return
                self.send(method(request))
            except Exception as error:
                self.send({"error": str(error)}, status=409)

        def log_message(self, *args):
            pass
    return Handler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=int)
    parser.add_argument("--kaon", action="store_true")
    parser.add_argument("--port", type=int, default=8767)
    args = parser.parse_args()
    if not 0 <= args.run <= 99999:
        parser.error("run must be 0–99999")
    if not BINARY.is_file():
        parser.error("build BH2_MIP_PREVIEW with CMake first")
    app = App(args.run, "K" if args.kaon else "Pi")
    server = ThreadingHTTPServer(("127.0.0.1", args.port), handler(app))
    print(f"run{args.run:05d} {app.particle} BH2 MIP / conda {sys.executable}", flush=True)
    print(f"open http://localhost:{args.port}/?run={args.run}  (Ctrl-C to stop)", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
