#!/usr/bin/env python3
"""BH2 ADC MIP and T0 tuning over a localhost SSH tunnel."""
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
from urllib.parse import parse_qs, urlparse

import numpy as np
import uproot

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from common.lib import config, fit_overrides  # noqa: E402
from common.lib.calib_state_io import locked_path, write_json  # noqa: E402

HERE = Path(__file__).resolve().parent
BINARY = PROJECT / "bin" / "BH2_MIP_PREVIEW"
T0_BINARY = PROJECT / "bin" / "T0_OFFSET_PREVIEW"
PLOTLY = Path(np.__file__).parents[1] / "plotly/package_data/plotly.min.js"


def read_json(path, default=None):
    if not path.exists():
        return copy.deepcopy(default)
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


class App:
    def __init__(self, run, particle="Pi", mode="hdprm"):
        if mode not in ("hdprm", "t0"):
            raise ValueError("unknown Hodo tuning mode")
        self.run, self.particle, self.mode = int(run), particle, mode
        self.expected_step = "hodo_t0" if mode == "t0" else "hodo_hdprm"
        self.stage_label = "T0 offset" if mode == "t0" else "HDPRM"
        self.binary = T0_BINARY if mode == "t0" else BINARY
        self.override_key = "gaussian_fit" if mode == "t0" else "mip_fit"
        self.validate_settings = fit_overrides.validate_t0_fit if mode == "t0" else fit_overrides.validate_mip_fit
        self.base = Path(config.OUTPUT_DIR) / "calib" / f"run{run:05d}"
        self.work = Path(config.WORK_DIR) / "mip_web" / f"run{run:05d}" / particle
        if mode == "t0":
            self.work = self.work / "t0"
        self.work.mkdir(parents=True, exist_ok=True)
        self.staging_path = self.base / f"{'t0' if mode == 't0' else 'mip'}_web_{particle}.json"
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
                   param / "USER" / config.SUB_DIR / f"UserParam_run{self.run:05d}", self.binary,
                   PROJECT / "common/include/params.h"]
        fit_root = self.fit_output()
        fit_stat = fit_root.stat() if fit_root.is_file() else None
        record = [str(resolved), stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns,
                  [(str(path), digest(path)) for path in sources],
                  [fit_stat.st_size, fit_stat.st_mtime_ns] if fit_stat else None]
        return hashlib.sha256(json.dumps(record).encode()).hexdigest()

    def channel(self, key):
        match = re.fullmatch(r"bh2-(\d+)-(t0)" if self.mode == "t0" else r"bh2-(\d+)-(u|d)", str(key))
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
                job["status"] = "done" if hold.get("step") == self.expected_step and hold.get("rc") == 0 else "failed"
                job["message"] = "再実行完了・Slackの新しいQAを確認してください" if job["status"] == "done" else "再実行が確認待ちで停止しました。Slackとログを確認してください"
            elif state.get("status") in ("failed", "stopped", "awaiting_final_approval", "approved") and hold.get("id") != job["hold_id"]:
                job["status"] = "done" if (state.get("status") in ("awaiting_final_approval", "approved")
                                                  and state.get("step") == self.expected_step
                                                  and state.get("last_rc", state.get("rc")) == 0) else "failed"
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
        lane_path = self.base / "lanes/hodo/calib_dispatch_state.json"
        paths = [lane_path] if lane_path.is_file() else [self.base / "calib_dispatch_state.json"]
        for path in paths:
            state = read_json(path, {})
            hold = state.get("hold") or {}
            held = state.get("status") == "held" and hold.get("step") == self.expected_step and not hold.get("decision")
            final_review = (state.get("status") == "awaiting_final_approval" and state.get("step") == self.expected_step
                            and state.get("finished") is True and state.get("rc") == 0
                            and type(state.get("index")) is int and state["index"] >= 0)
            if held or final_review:
                import yaml
                session_path = path.parent / "session.yml" if path == lane_path else self.base / "calib_session.yml"
                if not session_path.is_file():
                    session_path = self.base / "calib_session.yml"
                session = yaml.safe_load(session_path.read_text()) or {}
                if session.get("auto") or state.get("approval_mode") == "auto":
                    raise ValueError("自動進行のセッションではWriteできません。承認待ちを有効にしてください")
                from common.lib import calib_approval
                if final_review and calib_approval.effective_status(state) != calib_approval.STATUS_AWAITING:
                    raise ValueError("承認済みの結果は最終承認待ちに戻してから修正してください")
                _, suffix = calib_approval.particle_suffix(state)
                if suffix != self.particle:
                    raise ValueError("表示中の粒子と校正セッションが異なります")
                if final_review:
                    signature = hashlib.sha256(json.dumps(state, sort_keys=True).encode()).hexdigest()
                    hold = {"id": "final_review:"+signature, "source": "final_review", "step": self.expected_step,
                            "index": state["index"], "rc": 0}
                return path, state, hold
        raise ValueError(f"Writeは{self.stage_label}（{self.expected_step}）の確認待ちから行ってください。現在は閲覧・調整・仮保存のみ可能です")

    def public(self, doc):
        try:
            _, _, hold = self.gate()
            gate = {"hold_id": hold["id"], "reason": (f"最終承認待ちです。Writeすると{self.stage_label}だけを再実行し、新しいQAの確認へ戻ります"
                                                      if hold.get("source") == "final_review" else "")}
        except (ValueError, OSError) as error:
            gate = {"hold_id": None, "reason": str(error)}
        return {"run": self.run, "particle": self.particle, "mode": self.mode, "fingerprint": doc["fingerprint"],
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
        command = [str(self.binary), str(self.fit_input().resolve()), self.particle, str(segment), side, str(output)]
        if settings is not None:
            settings_path.write_text((fit_overrides.t0_fit_line if self.mode == "t0" else fit_overrides.mip_fit_line)(key, self.particle, settings) + "\n")
            # The preview uses the same particle-qualified key as the production reader.
            if curve_only and settings["function"] != "mean":
                command = [str(self.binary), "--curve", str(settings_path), f"{key}:{self.particle}", str(output)]
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

    def fit_output(self):
        kind = "T0_Offset" if self.mode == "t0" else "BH2_HDPRM"
        return Path(config.OUTPUT_DIR) / "root" / f"run{self.run:05d}" / f"run{self.run:05d}_{kind}_{self.particle}.root"

    def histogram_name(self, segment, side):
        return f"T0_seg{segment}_TimeOffset_{self.particle}" if self.mode == "t0" else f"BH2_ADC_seg{segment}{side}_{self.particle}"

    def channel_keys(self):
        return [f"bh2-{segment}-{side}" for segment in range(config.NUM_OF_CH["BH2"])
                for side in (("t0",) if self.mode == "t0" else ("u", "d"))]

    def fit_input(self):
        output = self.fit_output()
        if not output.is_file():
            raise ValueError(f"{self.stage_label}のFit結果がありません。先に当該段階を実行してください")
        with uproot.open(output) as source:
            required = [self.histogram_name(*self.channel(key)) for key in self.channel_keys()]
            if any(name not in source for name in required):
                raise ValueError(f"この{self.stage_label}結果には元の{'T0' if self.mode == 't0' else 'ADC'}ヒストグラムが保存されていません。当該段階を再実行してから開いてください")
        return output

    def overview(self, request):
        with self.lock, locked_path(self.staging_path):
            doc = self.document()
            self.require_current(request, doc)
            offset = int(request.get("page", 0)) * 12
            keys = self.channel_keys()
            if offset < 0 or offset >= len(keys):
                raise ValueError("page is invalid")
            with uproot.open(self.fit_input()) as source:
                channels = []
                for key in keys[offset:offset+12]:
                    segment, side = self.channel(key)
                    name = self.histogram_name(segment, side)
                    counts, edges = source[name].to_numpy()
                    if not np.all(np.isfinite(counts)):
                        raise ValueError(f"nonfinite histogram: {name}")
                    rebin = 3 if self.mode == "t0" else 2
                    used = len(counts)//rebin*rebin
                    counts = counts[:used].reshape(-1, rebin).sum(axis=1)
                    edges = edges[:used+1:rebin]
                    channels.append({"key": key, "label": f"BH2 seg{segment} {side}",
                                     "x": ((edges[:-1] + edges[1:])/2).tolist(), "counts": counts.tolist(),
                                     "bounds": [float(edges[0]), float(edges[-1])],
                                     "automatic": self.automatic(key, doc["fingerprint"])})
            self.require_current(request, self.document())
            return {**self.public(doc), "channels": channels, "pages": (len(keys)+11)//12}

    def invalid_snapshot(self, key, parameters, bounds, function):
        """Keep a bad saved fit editable, without drawing or adopting an invented fit."""
        segment, side = self.channel(key)
        with uproot.open(self.fit_input()) as source:
            counts, edges = source[self.histogram_name(segment, side)].to_numpy()
        rebin = 3 if self.mode == "t0" else 2
        upper = float(edges[len(counts)//rebin*rebin])
        lower = float(edges[0])
        span = upper-lower
        selected = list(bounds) if len(bounds) == 2 and np.all(np.isfinite(bounds)) and bounds[0] < bounds[1] else [lower, upper]
        selected = [max(lower, selected[0]), min(upper, selected[1])]
        if selected[0] >= selected[1]:
            selected = [lower, upper]
        centers = .5*(edges[:-1]+edges[1:])
        peak = float(centers[int(np.argmax(counts))])
        size = 3 if self.mode == "t0" else 4
        seeds = [max(1., float(np.max(counts))*span*.1), peak, max(float(edges[1]-edges[0]), span*.05), 0.][:size]
        for index, value in enumerate(parameters[:size]):
            if np.isfinite(value) and (index not in (0, 2) or value > 0.) and (index != 3 or value >= 0.):
                seeds[index] = float(value)
        seeds[1] = min(selected[1], max(selected[0], seeds[1]))
        return {"parameters": seeds, "range": selected, "function": function if function in ("gaussian", "landau", "mean") else "gaussian",
                "status": -2, "covariance": -1, "ndf": 0, "success": False, "curve": [], "unusable": True,
                "origin": "前回Fitの保存値が不正です。曲線は表示せず、調整用の初期値だけを用意しました。Fitまたは集計を実行してください"}

    def automatic(self, key, fingerprint):
        cache_key = (fingerprint, key)
        if cache_key not in self.cache:
            fit_root = self.fit_output()
            snapshot = None
            if fit_root.is_file():
                with uproot.open(fit_root) as output:
                    tree_name = "t0_fit" if self.mode == "t0" else "mip_fit"
                    if tree_name in output:
                        segment, side = self.channel(key)
                        columns = output[tree_name].arrays(library="np")
                        for index, row_seg in enumerate(columns["seg"]):
                            if row_seg == segment and (self.mode == "t0" or columns["side"][index] == (0 if side == "U" else 1)):
                                parameters = columns["parameters"][index].tolist()
                                bounds = columns["range"][index].tolist()
                                function = "gaussian" if self.mode == "t0" else str(columns["function"][index])
                                usable = (len(parameters) == (3 if self.mode == "t0" else 4)
                                          and np.all(np.isfinite(parameters+bounds)) and bounds[0] < bounds[1]
                                          and parameters[0] > 0 and parameters[2] > 0
                                          and function in ("gaussian", "landau", "mean")
                                          and (len(parameters) == 3 or parameters[3] >= 0.))
                                if usable:
                                    snapshot = {"parameters": parameters, "range": bounds, "function": function,
                                                "success": bool(columns["success"][index]), "status": int(columns["status"][index]),
                                                "covariance": int(columns["covariance"][index]), "ndf": int(columns["ndf"][index]),
                                                "origin": "前回実行の統計値（平均）" if function == "mean" else "前回実行のFit結果"}
                                    if "statistics" in columns and bool(columns["statistics_valid"][index]):
                                        values = columns["statistics"][index].tolist()
                                        snapshot["statistics"] = {"valid": True, "window": values[:2], "range": values[2:4],
                                            **dict(zip(("mean", "median", "window_median", "entries", "selected_entries", "error"), values[4:]))}
                                        if function == "mean":
                                            snapshot["range"] = snapshot["statistics"]["range"]
                                    elif snapshot["status"] < 0:
                                        snapshot["origin"] = ("Fit未実施（低統計）。初期値を調整してFitしてください" if self.mode == "t0"
                                                              else "Fit未実施（低統計）。詳細な統計値の表示にはHDPRMの再実行が必要です")
                                    curve_settings = {"parameters": parameters, "range": bounds,
                                                      "function": function, "range_mode": "fixed"}
                                    snapshot["curve"] = (self.execute_preview(key, curve_settings, curve_only=True)["curve"]
                                                         if snapshot["status"] >= 0 and function != "mean" else [])
                                else:
                                    snapshot = self.invalid_snapshot(key, parameters, bounds, function)
                                break
            if snapshot is None:
                raise ValueError(f"保存済みの{self.stage_label} Fit結果を読み出せません。当該段階を再実行してから開いてください")
            self.cache[cache_key] = snapshot
        return copy.deepcopy(self.cache[cache_key])

    def edit(self, request):
        with self.lock, locked_path(self.staging_path):
            doc = self.document()
            self.require_current(request, doc)
            key = request["key"]
            self.channel(key)
            result = self.automatic(key, doc["fingerprint"])
            entry = fit_overrides.load(self.run).get("fit", {}).get("BH2", {}).get(self.mode, {}).get(key, {})
            committed = entry.get(self.override_key, {}).get(self.particle)
            automatic = self.initial_settings(result)
            initial = doc["staged"].get(key, {}).get("settings") or committed or automatic
            self.require_current(request, self.document())
            return {**self.public(doc), "automatic": result, "initial": initial,
                    "baseline": committed or automatic, "automatic_note": result["origin"]}

    def initial_settings(self, result):
        return {"parameters": result["parameters"], "range": result["range"],
                "function": "gaussian" if self.mode == "t0" else "mean" if result["function"] == "mean" else "auto", "range_mode": "auto"}

    def fit(self, request, curve_only=False):
        with self.lock, locked_path(self.staging_path):
            doc = self.document()
            self.require_current(request, doc)
            key = request["key"]
            settings = self.validate_settings(request["settings"])
            if self.mode == "t0":
                segment, side = self.channel(key)
                with uproot.open(self.fit_input()) as source:
                    edges = source[self.histogram_name(segment, side)].axis().edges()
                    used = (len(edges)-1)//3*3
                    bounds = [float(edges[0]), float(edges[used])]
                if not bounds[0] <= settings["range"][0] < settings["range"][1] <= bounds[1]:
                    raise ValueError(f"T0のFit範囲は{bounds[0]}〜{bounds[1]} ns内で指定してください")
            elif not 0 <= settings["range"][0] < settings["range"][1] <= 4096:
                raise ValueError("Fit範囲はADCヒストグラム内（0–4096）で指定してください")
            if curve_only:
                effective = dict(settings)
                if effective["function"] == "auto":
                    effective["function"] = request.get("auto_function", "gaussian")
                result = self.execute_preview(key, effective, curve_only=True)
                self.require_current(request, self.document())
                return result
            result = self.execute_preview(key, settings)
            if settings["function"] == "mean" and not (result.get("statistics") or {}).get("valid"):
                result["success"] = False
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
                    raise ValueError("成功したFitまたは有効な集計を実行してから仮保存してください")
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
                entry = fit_overrides.load(self.run).get("fit", {}).get("BH2", {}).get(self.mode, {}).get(key, {})
                before = entry.get(self.override_key, {}).get(self.particle)
                if before is None:
                    automatic = self.automatic(key, doc["fingerprint"])
                    before = self.initial_settings(automatic)
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
                rows.append({"key": key, "changes": changes, "statistics": row["result"].get("statistics")})
            return {**self.public(doc), "rows": rows}

    def reopen_final_review(self, path, state, context):
        """An explicit Write reopens this unapproved final stage for refitting."""
        if context.get("source") != "final_review":
            return context
        sys.path.insert(0, str(PROJECT / "calibration/scripts"))
        import calib_dispatch
        hold = calib_dispatch.make_hold("web_refit_review", f"Web Write: refit {self.stage_label} before final approval",
                                        self.expected_step, state["index"], 0)
        lane = path.parent.name if path.parent.parent.name == "lanes" else ""
        if lane:
            hold["lane"] = lane
        reopened = copy.deepcopy(state)
        reopened.update(status="held", finished=False, hold=hold, gate_waiter=False,
                        hold_requested=False, stop_requested=False, updated_at=calib_dispatch.utc_now())
        reopened.pop("gate_waiter_pid", None)
        write_json(path, reopened)
        calib_dispatch.append_hold_journal(self.run, {"event": "web_refit_review", "at": hold["at"],
            "hold_id": hold["id"], "step": self.expected_step, "index": state["index"], "lane": lane,
            "source": hold["source"], "reason": hold["reason"], "final_review_id": context["id"]})
        calib_dispatch.append_progress(self.run, "web_refit_review", step=self.expected_step,
                                       index=state["index"], status="held", hold_source="web_refit_review")
        return hold

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
                channels = overrides.setdefault("fit", {}).setdefault("BH2", {}).setdefault(self.mode, {})
                for key, row in doc["staged"].items():
                    channels.setdefault(key, {}).setdefault(self.override_key, {})[self.particle] = self.validate_settings(row["settings"])
                target = fit_overrides.overrides_path(self.run)
                if target.exists():
                    backup = self.work / (time.strftime("overrides_%Y%m%dT%H%M%S_") + secrets.token_hex(4) + ".json")
                    shutil.copy2(target, backup)
                fit_overrides.save(self.run, overrides)
                doc["job"] = {"status": "running", "hold_id": hold["id"], "state_path": str(state_path),
                              "started_at": time.time(), "message": f"{self.stage_label}を再実行中。完了後にSlackの新しいQAで確認待ちになります"}
                write_json(self.staging_path, doc)
                # Refit uses the existing lane's decode mode and always returns to its gate.
                from common.lib import calib_ops
                try:
                    hold = self.reopen_final_review(state_path, state, hold)
                    doc["job"]["hold_id"] = hold["id"]
                    write_json(self.staging_path, doc)
                    accepted, message = calib_ops.run_dispatch("retry", self.run, hold_id=hold["id"], background_resume=True)
                    if not accepted:
                        raise RuntimeError(message)
                except BaseException as error:
                    doc["job"].update(status="failed", message=f"設定は保存済み。再実行の開始に失敗: {error}")
                    doc["fingerprint"] = self.fingerprint()
                    write_json(self.staging_path, doc)
                return self.public(doc)



class Portal:
    """Keep independent run/particle/stage workspaces behind one localhost endpoint."""
    def __init__(self, default_run=None, default_particle=None):
        self.default_run = default_run
        self.default_particle = default_particle
        self.applications = {}
        self.result_headers = {}
        self.lock = threading.Lock()
        self.preview_lock = threading.Lock()

    def available(self):
        rows = []
        folder = Path(config.OUTPUT_DIR) / "root"
        pattern = re.compile(r"run(\d{5})_(BH2_HDPRM|T0_Offset)_(Pi|K)\.root")
        for path in folder.glob("run*/run*.root"):
            match = pattern.fullmatch(path.name)
            if not match or path.parent.name != "run"+match[1]:
                continue
            run = int(match[1])
            if not (Path(config.DATA_DIR) / f"run{run:05d}_Hodo.root").is_file():
                continue
            mode = "t0" if match[2] == "T0_Offset" else "hdprm"
            try:
                stat = path.stat()
                stamp = (str(path), stat.st_mtime_ns, stat.st_size)
                with self.lock:
                    readiness = self.result_headers.get(stamp)
                if readiness is None:
                    particle = match[3]
                    required = ([f"T0_seg{segment}_TimeOffset_{particle}" for segment in range(config.NUM_OF_CH["BH2"])]
                                if mode == "t0" else [f"BH2_ADC_seg{segment}{side}_{particle}"
                                    for segment in range(config.NUM_OF_CH["BH2"]) for side in ("U", "D")])
                    with uproot.open(path) as source:
                        ready = all(name in source for name in required)
                    readiness = (ready, "" if ready else "元ヒストグラムなし・当該段階の再実行が必要")
                    with self.lock:
                        self.result_headers[stamp] = readiness
            except (OSError, ValueError):
                readiness = (False, "結果ファイルを読み取れません")
            rows.append({"run": run, "particle": match[3], "mode": mode, "ready": readiness[0], "reason": readiness[1]})
        return sorted(rows, key=lambda row: (-row["run"], row["particle"], row["mode"]))

    def infer_particle(self, run, mode):
        base = Path(config.OUTPUT_DIR) / "calib" / f"run{run:05d}"
        lane = base / "lanes/hodo/calib_dispatch_state.json"
        state = read_json(lane if lane.is_file() else base / "calib_dispatch_state.json", {})
        plan = state.get("plan") or {}
        if state.get("calib_particle") in ("pion", "kaon") or plan.get("calib_particle") in ("pion", "kaon") or state.get("kaon") is True or plan.get("kaon") is True:
            from common.lib.calib_approval import particle_suffix
            return particle_suffix(state)[1]
        kind = "T0_Offset" if mode == "t0" else "BH2_HDPRM"
        choices = [particle for particle in ("Pi", "K") if
                   (Path(config.OUTPUT_DIR) / "root" / f"run{run:05d}" / f"run{run:05d}_{kind}_{particle}.root").is_file()]
        if len(choices) == 1:
            return choices[0]
        if self.default_particle and run == self.default_run:
            return self.default_particle
        if len(choices) > 1:
            raise ValueError("PiとKの結果があります。入口の一覧から粒子を選んでください")
        raise ValueError("この段階の保存済みFit結果がありません。先に較正を実行してください")

    def resolve(self, query):
        if any(len(query.get(key, [])) > 1 for key in ("run", "mode", "particle")):
            raise ValueError("対象の指定が重複しています")
        value = query.get("run", [str(self.default_run) if self.default_run is not None else ""])[0]
        if not re.fullmatch(r"\d{1,5}", value):
            raise ValueError("対象runを一覧またはSlackのリンクから選んでください")
        run = int(value)
        mode = query.get("mode", ["hdprm"])[0]
        if mode not in ("hdprm", "t0"):
            raise ValueError("unknown Hodo tuning mode")
        if not (Path(config.DATA_DIR) / f"run{run:05d}_Hodo.root").is_file():
            raise ValueError(f"run{run:05d}のHodo decodeがありません")
        particle = query.get("particle", [None])[0] or self.infer_particle(run, mode)
        if particle not in ("Pi", "K"):
            raise ValueError("粒子はPiまたはKを選んでください")
        key = (run, particle, mode)
        with self.lock:
            if key not in self.applications:
                app = App(run, particle, mode)
                app.preview_lock = self.preview_lock
                self.applications[key] = app
            return self.applications[key]


def handler(app=None, t0_app=None, portal=None):
    applications = {"hdprm": app}
    if t0_app is not None:
        applications["t0"] = t0_app
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

        def application(self):
            query = parse_qs(urlparse(self.path).query)
            if portal is not None:
                return portal.resolve(query)
            mode = query.get("mode", ["hdprm"])[0]
            if mode not in applications:
                raise ValueError("unknown Hodo tuning mode")
            return applications[mode]

        def do_GET(self):
            if not self.valid_host():
                return
            try:
                path = urlparse(self.path).path
                query = parse_qs(urlparse(self.path).query)
                if portal is not None and path == "/portal-health":
                    self.send({"service": "e72_calib_portal", "version": 1})
                    return
                if portal is not None and path == "/runs":
                    self.send({"runs": portal.available()})
                    return
                if portal is not None and path == "/" and "run" not in query and portal.default_run is None:
                    self.send((HERE / "calib_portal.html").read_bytes(), "text/html; charset=utf-8")
                    return
                if path == "/plotly.min.js":
                    self.send(PLOTLY.read_bytes(), "application/javascript")
                    return
                try:
                    app = self.application()
                except (ValueError, OSError) as error:
                    if portal is not None and path == "/":
                        import html
                        page = (HERE / "calib_portal.html").read_text().replace("__ERROR__", html.escape(str(error)), 1)
                        self.send(page.encode(), "text/html; charset=utf-8")
                        return
                    raise
                if urlparse(self.path).path == "/state":
                    self.send(app.state())
                elif urlparse(self.path).path == "/":
                    page = (HERE / "mip_web.html").read_text().replace("__TOKEN__", app.token).replace("__RUN__", str(app.run)).replace("__PARTICLE__", app.particle)
                    self.send(page.encode(), "text/html; charset=utf-8")
                else:
                    self.send({"error": "not found"}, status=404)
            except Exception as error:
                self.send({"error": str(error)}, status=409)

        def do_POST(self):
            if not self.valid_host():
                return
            try:
                app = self.application()
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
    parser.add_argument("run", type=int, nargs="?")
    parser.add_argument("--kaon", action="store_true")
    parser.add_argument("--port", type=int, default=8767)
    args = parser.parse_args()
    if args.run is not None and not 0 <= args.run <= 99999:
        parser.error("run must be 0–99999")
    if not BINARY.is_file():
        parser.error("build BH2_MIP_PREVIEW with CMake first")
    if not T0_BINARY.is_file():
        parser.error("build T0_OFFSET_PREVIEW with CMake first")
    portal = Portal(args.run, "K" if args.kaon else None)
    server = ThreadingHTTPServer(("127.0.0.1", args.port), handler(portal=portal))
    query = f"?run={args.run}" if args.run is not None else ""
    if args.kaon:
        query += ("&" if query else "?") + "particle=K"
    print(f"Calibration review / conda {sys.executable}", flush=True)
    print(f"open http://localhost:{args.port}/{query}  (Ctrl-C to stop)", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
