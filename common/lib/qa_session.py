"""Bind QA PDFs to one fit execution via content-hash snapshots.

begin() records expected PDF hashes before a fit. verify_sources() requires
each file to exist and differ from that snapshot (or be newly created).
mtime alone is never used for pass/fail.

Successful outputs are sealed under qa_manifests/<exec_id>/ so resend can
re-read an immutable copy and refuse hash mismatch.
"""
from __future__ import annotations

import hashlib
import json
import secrets
import shutil
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
def seal_dir(session: dict[str, Any]) -> Path:
    """Directory for immutable sealed PDFs of this exec_id."""
    return manifest_dir(int(session["run"])) / str(session["exec_id"])


#_____________________________________________________________________________
def seal_outputs(session: dict[str, Any], output_pdfs) -> list[Path]:
    """Copy output PDFs under qa_manifests/<exec_id>/; return sealed paths."""
    dest_dir = seal_dir(session)
    dest_dir.mkdir(parents=True, exist_ok=True)
    sealed: list[Path] = []
    for path in output_pdfs:
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(f"cannot seal missing PDF: {path}")
        dest = dest_dir / path.name
        shutil.copy2(path, dest)
        sealed.append(dest)
    return sealed


#_____________________________________________________________________________
def attach_params(session: dict[str, Any], paths, *, role: str = "fit_input") -> list[dict[str, Any]]:
    """Snapshot parameter files onto the session (sha at call time)."""
    entries: list[dict[str, Any]] = []
    for path in paths:
        path = Path(path)
        entries.append({
            "name": path.name,
            "path": str(path),
            "sha": file_sha256(path) if path.is_file() else None,
            "role": role,
        })
    if role == "fit_input":
        session["params_fit_input"] = entries
    elif role == "written":
        session["params_written"] = entries
    else:
        session[f"params_{role}"] = entries
    return entries


#_____________________________________________________________________________
def write_manifest(
    session: dict[str, Any],
    *,
    source_hashes: dict[str, str],
    output_pdfs,
    param_paths=None,
    params_written=None,
) -> Path:
    """Write qa_manifests/<exec_id>.json; return the path.

    output_pdfs should be sealed paths when available. params come from
    session['params_fit_input'] (preferred) or param_paths (role fit_input).
    """
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
    params = list(session.get("params_fit_input") or [])
    if not params and param_paths:
        for path in param_paths:
            path = Path(path)
            params.append({
                "name": path.name,
                "path": str(path),
                "sha": file_sha256(path) if path.is_file() else None,
                "role": "fit_input",
            })
    written = list(session.get("params_written") or [])
    if params_written:
        written = []
        for path in params_written:
            path = Path(path)
            written.append({
                "name": path.name,
                "path": str(path),
                "sha": file_sha256(path) if path.is_file() else None,
                "role": "written",
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
        "params_written": written,
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
def resolve_sealed_outputs(
    manifest: dict[str, Any],
    *,
    run: int,
    mode: str | None = None,
    particle: str | None = None,
) -> tuple[bool, str, list[Path]]:
    """Validate run/mode/particle and on-disk sha against manifest outputs."""
    if int(manifest.get("run", -1)) != int(run):
        return (
            False,
            f"manifest run mismatch: manifest={manifest.get('run')} cli={run}",
            [],
        )
    meta = manifest.get("meta") or {}
    if mode is not None and meta.get("mode") is not None:
        if str(meta.get("mode")) != str(mode):
            return (
                False,
                f"manifest mode mismatch: manifest={meta.get('mode')} cli={mode}",
                [],
            )
    if particle is not None and meta.get("particle") is not None:
        if str(meta.get("particle")) != str(particle):
            return (
                False,
                (
                    f"manifest particle mismatch: "
                    f"manifest={meta.get('particle')} cli={particle}"
                ),
                [],
            )
    outputs = manifest.get("outputs") or []
    if not outputs:
        return False, "manifest has no outputs", []
    paths: list[Path] = []
    for entry in outputs:
        path = Path(entry.get("path") or "")
        expected_sha = entry.get("sha")
        if not path.is_file():
            return False, f"sealed PDF missing: {path.name or path}", []
        actual = file_sha256(path)
        if expected_sha is not None and actual != expected_sha:
            return (
                False,
                f"sealed PDF sha mismatch: {path.name} "
                f"(manifest={expected_sha} disk={actual})",
                [],
            )
        paths.append(path)
    return True, "ok", paths


#_____________________________________________________________________________
def resolve_sealed_output(
    manifest: dict[str, Any],
    *,
    run: int,
    mode: str | None = None,
    particle: str | None = None,
) -> tuple[bool, str, Path | None]:
    """Single-PDF wrapper around resolve_sealed_outputs."""
    ok, message, paths = resolve_sealed_outputs(
        manifest, run=run, mode=mode, particle=particle,
    )
    if not ok:
        return False, message, None
    return True, message, paths[0]


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
