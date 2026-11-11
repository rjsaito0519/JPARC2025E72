"""Shared isolated tree for calib_dispatch unit tests (#39).

Copies dispatcher deps and stubs site/Slack/ROOT-facing helpers so an empty
PYTHONPATH tree can import and run calib_dispatch without the real repo.
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any


#_____________________________________________________________________________
def populate_dispatch_lib(
    repo: Path,
    root: Path,
    *,
    analyzer_dir: Path | None = None,
    extra_config_lines: str = "",
) -> dict[str, Any]:
    """Create common/lib + termcolor for an isolated calib_dispatch tree.

    Returns paths: scripts, lib, data, decode, scratch, out, param.
    """
    root = Path(root)
    repo = Path(repo)
    scripts = root / "calibration" / "scripts"
    scripts.mkdir(parents=True, exist_ok=True)
    lib = root / "common" / "lib"
    lib.mkdir(parents=True, exist_ok=True)
    (lib / "__init__.py").touch()

    data = root / "data"
    decode = root / "decode"
    scratch = root / "scratch"
    out = root / "out"
    param = root / "param"
    for path in (data, decode, scratch, out, param):
        path.mkdir(parents=True, exist_ok=True)

    config_bits = [
        "from pathlib import Path\n",
        f"OUTPUT_DIR = Path({str(out)!r})\n",
        f"PARAM_DIR = Path({str(param)!r})\n",
        f"DATA_DIR = Path({str(data)!r})\n",
        f"DECODE_DIR = Path({str(decode)!r})\n",
        f"SCRATCH_DIR = Path({str(scratch)!r})\n",
        "SUB_DIR = 'e72'\n",
    ]
    if analyzer_dir is not None:
        config_bits.insert(
            3, f"ANALYZER_DIR = Path({str(analyzer_dir)!r})\n",
        )
    if extra_config_lines:
        config_bits.append(extra_config_lines if extra_config_lines.endswith("\n")
                           else extra_config_lines + "\n")
    (lib / "config.py").write_text("".join(config_bits), encoding="utf-8")

    (lib / "input_check.py").write_text(
        "def check_input(*a, **k): return True\n"
        "def stale_params(*a, **k): return []\n",
        encoding="utf-8",
    )
    (lib / "param.py").write_text("def backup_file(*a): pass\n", encoding="utf-8")
    (root / "termcolor.py").write_text(
        "def colored(text, *a, **kw): return text\n", encoding="utf-8",
    )
    (lib / "channel_flags.py").write_text(
        "def open_needs_review(*a, **k):\n"
        "    return []\n"
        "def write_blocked_by_flags(*a, **k):\n"
        "    return False, ''\n",
        encoding="utf-8",
    )
    (lib / "final_qa_index.py").write_text(
        "def index_summary(*a, **k):\n"
        "    return {'n_pdf': 0, 'pdfs': []}\n"
        "def index_from_subject(run, subject):\n"
        "    return {'run': int(run), 'n_pdf': 0, 'pdfs': [], "
        "'source': 'subject_sealed'}\n",
        encoding="utf-8",
    )
    (lib / "slack_notify.py").write_text(
        "def try_notify_failure(*a, **k):\n"
        "    return True, 'slack skip: disabled'\n",
        encoding="utf-8",
    )
    (lib / "slack_bot_ensure.py").write_text(
        "TMUX_SESSION = 'e72_calib_slack_bot'\n"
        "def ensure_slack_bot(*, required=False):\n"
        "    return True, 'slack bot ensure stub (tests)'\n",
        encoding="utf-8",
    )

    shutil.copy(repo / "calibration" / "scripts" / "calib_dispatch.py", scripts)
    for name in (
        "runner.py", "calib_approval.py", "qa_session.py",
        "calib_anomaly.py", "calib_progress.py",
        "calib_lanes.py", "fit_validity.py", "calib_list.py",
    ):
        shutil.copy(repo / "common" / "lib" / name, lib)
    # Minimal stubs for ./calib.py imports used by isolated CLI tests.
    if not (lib / "hodo_btof_monitor.py").is_file():
        (lib / "hodo_btof_monitor.py").write_text(
            "from pathlib import Path\n"
            "from . import config\n"
            "def hodo_root(run):\n"
            "    return Path(config.DATA_DIR) / f'run{int(run):05d}_Hodo.root'\n",
            encoding="utf-8",
        )
    if not (lib / "calib_runlist.py").is_file():
        shutil.copy(repo / "common" / "lib" / "calib_runlist.py", lib)
    if not (lib / "decode_run.py").is_file():
        shutil.copy(repo / "common" / "lib" / "decode_run.py", lib)

    return {
        "scripts": scripts,
        "lib": lib,
        "data": data,
        "decode": decode,
        "scratch": scratch,
        "out": out,
        "param": param,
    }
