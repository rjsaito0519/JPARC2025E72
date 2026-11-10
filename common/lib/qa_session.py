"""Bind QA PDFs to one fit execution via content-hash snapshots.

begin() records expected PDF hashes before a fit. verify_sources() requires
each file to exist and differ from that snapshot (or be newly created).
mtime alone is never used for pass/fail.
"""
from __future__ import annotations

import hashlib
import json
import secrets
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from . import config


#_____________________________________________________________________________
def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def sha256_hex(data: bytes | str, n: int = 16) -> str:
    if isinstance(data, str):
        data = data.encode()
    return hashlib.sha256(data).hexdigest()[:n]


#_____________________________________________________________________________
def file_sha256(path: Path | str, n: int = 16) -> str:
    return sha256_hex(Path(path).read_bytes(), n=n)


#_____________________________________________________________________________
def input_token(path: Path | str) -> str:
    """Id of an input ROOT (name + size + mtime_ns)."""
    path = Path(path)
    st = path.stat()
    return f"{path.name}:{st.st_size}:{st.st_mtime_ns}"


#_____________________________________________________________________________
def snapshot_hash(path: Path | str) -> Optional[str]:
    path = Path(path)
    if not path.is_file():
        return None
    return file_sha256(path)


#_____________________________________________________________________________
def begin(
    run: int,
    label: str,
    expected_paths,
    *,
    input_path: Path | str | None = None,
    meta: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Start a QA session: snapshot hashes of expected PDFs before the fit."""
    expected = [Path(p) for p in expected_paths]
    started = _utc_now()
    exec_id = f"{started.replace(':', '').replace('-', '')}_{secrets.token_hex(4)}"
    tok = None
    if input_path is not None:
        ip = Path(input_path)
        if ip.is_file():
            tok = input_token(ip)
    return {
        "exec_id": exec_id,
        "run": int(run),
        "label": str(label),
        "meta": dict(meta or {}),
        "input_token": tok,
        "expected": [str(p) for p in expected],
        "before": {p.name: snapshot_hash(p) for p in expected},
        "started_at": started,
    }


#_____________________________________________________________________________
def verify_sources(session: dict[str, Any]) -> tuple[bool, str, dict[str, str]]:
    """Require every expected PDF to exist and differ from the begin() snapshot."""
    hashes: dict[str, str] = {}
    missing: list[str] = []
    stale: list[str] = []
    before = session.get("before") or {}
    for path_str in session.get("expected") or []:
        path = Path(path_str)
        name = path.name
        if not path.is_file():
            missing.append(name)
            continue
        digest = file_sha256(path)
        hashes[name] = digest
        prev = before.get(name)
        if prev is not None and prev == digest:
            stale.append(name)
    if missing:
        return False, f"missing PDF(s): {', '.join(missing)}", hashes
    if stale:
        return (
            False,
            f"stale PDF(s) unchanged since fit start: {', '.join(stale)}",
            hashes,
        )
    return True, "ok", hashes


#_____________________________________________________________________________
def manifest_dir(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / f"run{run:05d}" / "qa_manifests"


#_____________________________________________________________________________
def write_manifest(
    session: dict[str, Any],
    *,
    source_hashes: dict[str, str],
    output_pdfs,
    param_paths=None,
) -> Path:
    """Write qa_manifests/<exec_id>.json; return the path."""
    run = int(session["run"])
    out_dir = manifest_dir(run)
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for path in output_pdfs:
        path = Path(path)
        outputs.append({
            "name": path.name,
            "path": str(path),
            "sha": file_sha256(path) if path.is_file() else None,
        })
    params = []
    for path in (param_paths or []):
        path = Path(path)
        params.append({
            "name": path.name,
            "path": str(path),
            "sha": file_sha256(path) if path.is_file() else None,
        })
    sources = [
        {"name": name, "sha": source_hashes.get(name)}
        for name in sorted(source_hashes)
    ]
    payload = {
        "exec_id": session["exec_id"],
        "run": run,
        "label": session.get("label"),
        "meta": session.get("meta") or {},
        "input_token": session.get("input_token"),
        "sources": sources,
        "outputs": outputs,
        "params": params,
        "started_at": session.get("started_at"),
        "finished_at": _utc_now(),
    }
    path = out_dir / f"{session['exec_id']}.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


#_____________________________________________________________________________
def load_manifest(path: Path | str) -> dict[str, Any]:
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or not data.get("exec_id"):
        raise ValueError(f"invalid qa manifest: {path}")
    return data


#_____________________________________________________________________________
def manifest_post_fields(manifest: dict[str, Any], pdf_path: Path | str) -> dict[str, Any]:
    """Fields to store on slack_thread.json posts[key]."""
    pdf_path = Path(pdf_path)
    pdf_sha = None
    for entry in manifest.get("outputs") or []:
        if entry.get("name") == pdf_path.name and entry.get("sha"):
            pdf_sha = entry["sha"]
            break
    if pdf_sha is None and pdf_path.is_file():
        pdf_sha = file_sha256(pdf_path)
    return {
        "exec_id": manifest.get("exec_id"),
        "input_token": manifest.get("input_token"),
        "pdf_sha": pdf_sha,
        "manifest": f"{manifest.get('exec_id')}.json",
    }
