#!/usr/bin/env python3
"""
TPC calibration auto-iterator.

Outer schedule: hit, hit, trk, ... (2:1)

hit (TPCHitBcOut):
  recon -> phase(--ofs-update) -> recon -> offset
        -> recon -> phase -> recon -> offset

trk (tracking):
  recon -> offset -> recon -> drift
        -> phase(--debug)   # eval only

recon = tmux split-pane myrun.py -> poll stat JSON + fresh ROOT -> add.sh/add_trk.sh

Usage (recommended inside tmux session `tpc_iter`):
  ./calibration/scripts/iterate_tpc_calib.py --max-iter 6
  ./calibration/scripts/iterate_tpc_calib.py --max-iter 1 --dry-run
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from dataclasses import dataclass
from typing import Optional

import yaml

try:
    from discordwebhook import Discord
except ImportError:  # pragma: no cover
    Discord = None  # type: ignore

# Same webhook as runmanager/myrun.py
DISCORD_WEBHOOK = (
    "https://discord.com/api/webhooks/1195721561018220544/"
    "g8IqwcoYDotpcmXTrEgB3o-Bnc_GXsAU_lZzAm7zBZtWi_pjLi3z0LAAtbXxgoI48Tfa"
)

SCRIPT_DIR = Path(__file__).resolve().parent
# myanalysis root (== /home/had/sryuta/JPARC2025E72 via symlink)
MYANALYSIS_ROOT = SCRIPT_DIR.parent.parent
RUN_TPC = SCRIPT_DIR / "run_tpc.py"

HIT_YML_REL = Path("runlist/dst_tpchit_bcout.yml")
TRK_YML_REL = Path("runlist/dst_tpctracking.yml")
POLL_SEC = 30
RECON_TIMEOUT_SEC = 7 * 24 * 3600  # 7 days
MTIME_SLACK_SEC = 2.0  # clock / filesystem slack for ROOT mtime vs t0


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


@dataclass
class StatSnapshot:
    failed: list[str]
    unfinished: list[str]
    done_count: int
    expected_count: int

    @property
    def all_done(self) -> bool:
        return (
            self.expected_count > 0
            and not self.failed
            and not self.unfinished
            and self.done_count == self.expected_count
        )


@dataclass
class MyrunHandle:
    tag: str
    token: str
    pane_id: str
    ec_file: Path
    foreground: bool = False
    exit_code: Optional[int] = None


class IterOrchestrator:
    def __init__(
        self,
        analyzer_dir: Path,
        max_iter: int,
        dry_run: bool,
        poll_sec: int = POLL_SEC,
        recon_timeout_sec: int = RECON_TIMEOUT_SEC,
    ) -> None:
        self.analyzer_dir = analyzer_dir.resolve()
        self.max_iter = max_iter
        self.dry_run = dry_run
        self.poll_sec = poll_sec
        self.recon_timeout_sec = recon_timeout_sec
        self.group_root = self.analyzer_dir / "group_root"
        self.stat_dir = self.analyzer_dir / "runmanager" / "stat"
        self.log_dir = MYANALYSIS_ROOT / "results" / "tpc_calib_iter"
        self.log_path = self.log_dir / "iter_log.txt"
        self.discord = None
        if not dry_run and Discord is not None:
            self.discord = Discord(url=DISCORD_WEBHOOK)

        self.hit_yml = self.analyzer_dir / HIT_YML_REL
        self.trk_yml = self.analyzer_dir / TRK_YML_REL
        self._validate_paths()

        self.hit_runs, self.hit_roots = self._parse_runlist(self.hit_yml)
        self.trk_runs, self.trk_roots = self._parse_runlist(self.trk_yml)
        if not self.hit_runs or not self.trk_runs:
            raise SystemExit("[Error] No active RUN entries in hit/trk yml")
        self.hit_run0 = self.hit_runs[0]
        self.trk_run0 = self.trk_runs[0]

    def _validate_paths(self) -> None:
        for p in (self.hit_yml, self.trk_yml, RUN_TPC, self.group_root):
            if not p.exists():
                raise SystemExit(f"[Error] Missing path: {p}")
        for name in ("add.sh", "add_trk.sh", "myrun"):
            pass
        if not (self.group_root / "add.sh").is_file():
            raise SystemExit(f"[Error] Missing {self.group_root / 'add.sh'}")
        if not (self.group_root / "add_trk.sh").is_file():
            raise SystemExit(f"[Error] Missing {self.group_root / 'add_trk.sh'}")
        myrun = self.analyzer_dir / "runmanager" / "myrun.py"
        if not myrun.is_file():
            raise SystemExit(f"[Error] Missing {myrun}")

    @staticmethod
    def _parse_runlist(yml_path: Path) -> tuple[list[int], list[Path]]:
        with yml_path.open(encoding="utf-8") as f:
            data = yaml.safe_load(f)
        runs: list[int] = []
        roots: list[Path] = []
        run_block = data.get("RUN") or {}
        for key, cfg in run_block.items():
            if cfg is None:
                continue
            try:
                run_num = int(key)
            except (TypeError, ValueError):
                continue
            root = cfg.get("root")
            if not root:
                continue
            runs.append(run_num)
            roots.append(Path(root))
        return runs, roots

    def log(self, msg: str) -> None:
        line = f"[{now()}] {msg}"
        print(line, flush=True)
        if self.dry_run:
            return
        self.log_dir.mkdir(parents=True, exist_ok=True)
        with self.log_path.open("a", encoding="utf-8") as f:
            f.write(line + "\n")

    def notify(self, content: str) -> None:
        self.log(f"Discord: {content}")
        if self.dry_run or self.discord is None:
            return
        try:
            self.discord.post(content=content[:1900])
        except Exception as exc:  # noqa: BLE001
            self.log(f"[WARN] Discord post failed: {exc}")

    def fail(self, msg: str) -> None:
        self.notify(f"[tpc_iter] STOP: {msg}")
        raise SystemExit(f"[Error] {msg}")

    @staticmethod
    def kind_for_outer(i: int) -> str:
        """1-based outer index -> hit, hit, trk, ..."""
        return "trk" if (i % 3 == 0) else "hit"

    def yml_for(self, kind: str) -> Path:
        return self.hit_yml if kind == "hit" else self.trk_yml

    def roots_for(self, kind: str) -> list[Path]:
        return self.hit_roots if kind == "hit" else self.trk_roots

    def run0_for(self, kind: str) -> int:
        return self.hit_run0 if kind == "hit" else self.trk_run0

    def hadd_root_for(self, kind: str) -> Path:
        name = (
            "htof_calib_hsoff.root"
            if kind == "hit"
            else "htof_calib_tracking_hsoff.root"
        )
        return self.group_root / name

    def add_script_for(self, kind: str) -> Path:
        return self.group_root / ("add.sh" if kind == "hit" else "add_trk.sh")

    def stat_json_for(self, yml: Path) -> Path:
        return self.stat_dir / f"{yml.stem}.json"

    @staticmethod
    def expected_run_keys(yml: Path) -> list[str]:
        with yml.open(encoding="utf-8") as f:
            data = yaml.safe_load(f)
        keys: list[str] = []
        for key, cfg in (data.get("RUN") or {}).items():
            if cfg is None:
                continue
            if not cfg.get("root"):
                continue
            keys.append(str(key))
        return keys

    def _read_stat_info(self, path: Path) -> Optional[dict]:
        if not path.is_file():
            return None
        raw = path.read_text(encoding="utf-8").strip()
        if not raw:
            return None
        try:
            info = json.loads(raw)
        except json.JSONDecodeError:
            return None
        return info if isinstance(info, dict) else None

    def _analyze_stat(self, info: dict, expected_keys: list[str]) -> StatSnapshot:
        failed: list[str] = []
        unfinished: list[str] = []
        done_count = 0
        for key in expected_keys:
            item = info.get(key)
            if not isinstance(item, dict):
                unfinished.append(f"{key}:missing")
                continue
            stat = str(item.get("stat", ""))
            if "FAIL" in stat.upper():
                failed.append(f"{key}:{stat}")
            elif stat == "DONE":
                done_count += 1
            else:
                unfinished.append(f"{key}:{stat or '?'}")
        return StatSnapshot(
            failed=failed,
            unfinished=unfinished,
            done_count=done_count,
            expected_count=len(expected_keys),
        )

    def check_output_roots_fresh(
        self, kind: str, t0: float
    ) -> tuple[list[str], list[str]]:
        """Return (missing, stale) ROOT paths relative to recon start time t0."""
        if self.dry_run:
            return [], []
        missing: list[str] = []
        stale: list[str] = []
        cutoff = t0 - MTIME_SLACK_SEC
        for path in self.roots_for(kind):
            if not path.is_file():
                missing.append(str(path))
                continue
            if path.stat().st_mtime < cutoff:
                stale.append(str(path))
        return missing, stale

    def _myrun_ec_ready(self, handle: MyrunHandle) -> bool:
        if handle.foreground:
            return handle.exit_code is not None
        return handle.ec_file.is_file()

    def _myrun_read_ec(self, handle: MyrunHandle) -> int:
        if handle.foreground:
            return handle.exit_code if handle.exit_code is not None else 1
        try:
            return int(handle.ec_file.read_text(encoding="utf-8").strip())
        except (OSError, ValueError):
            return 1

    def wait_for_recon_complete(
        self,
        yml: Path,
        kind: str,
        t0: float,
        handle: Optional[MyrunHandle],
        tag: str,
    ) -> None:
        """Poll stat JSON until expected runs DONE and ROOT outputs are fresh."""
        if self.dry_run:
            self.log(
                f"[dry-run] would poll {self.stat_json_for(yml)} "
                f"every {self.poll_sec}s"
            )
            return

        stat_path = self.stat_json_for(yml)
        expected = self.expected_run_keys(yml)
        if not expected:
            self.fail(f"No RUN keys in {yml}")

        deadline = (
            time.time() + self.recon_timeout_sec
            if self.recon_timeout_sec > 0
            else None
        )
        last_snap: Optional[StatSnapshot] = None
        poll_n = 0
        myrun_done_at: Optional[float] = None
        myrun_grace_sec = max(3 * self.poll_sec, 90)

        self.log(
            f"Polling recon ({tag}): stat={stat_path.name} "
            f"runs={expected} interval={self.poll_sec}s"
        )

        while True:
            poll_n += 1
            info = self._read_stat_info(stat_path)
            if info is not None:
                snap = self._analyze_stat(info, expected)
                last_snap = snap
                if snap.failed:
                    self.fail(
                        f"DST FAILED in {stat_path.name} ({tag}): "
                        f"{', '.join(snap.failed)}"
                    )

            if handle is not None and self._myrun_ec_ready(handle):
                ec = self._myrun_read_ec(handle)
                if ec != 0:
                    snap_msg = ""
                    if last_snap is not None:
                        snap_msg = (
                            f"; stat done={last_snap.done_count}/"
                            f"{last_snap.expected_count}"
                        )
                    self.fail(f"myrun exit {ec} ({tag}){snap_msg}")
                if myrun_done_at is None:
                    myrun_done_at = time.time()

            if last_snap is not None and last_snap.all_done:
                missing, stale = self.check_output_roots_fresh(kind, t0)
                if not missing and not stale:
                    if handle is None or self._myrun_ec_ready(handle):
                        break
                    if poll_n % 10 == 0:
                        self.log(
                            f"recon ({tag}): stat/ROOT OK; "
                            "waiting for myrun process exit"
                        )
                elif poll_n == 1 or poll_n % 10 == 0:
                    kind_msg = "missing" if missing else "stale"
                    n_bad = len(missing) if missing else len(stale)
                    self.log(
                        f"recon ({tag}): stat DONE but {kind_msg} ROOT "
                        f"({n_bad}/{len(expected)})"
                    )
            elif poll_n == 1:
                self.log(f"recon ({tag}): stat not ready yet ({stat_path.name})")

            if myrun_done_at is not None:
                elapsed = time.time() - myrun_done_at
                if elapsed > myrun_grace_sec:
                    if last_snap is None or not last_snap.all_done:
                        self.fail(
                            f"myrun exit 0 but stat incomplete ({tag}): "
                            f"{last_snap.unfinished if last_snap else 'no stat'}"
                        )
                    missing, stale = self.check_output_roots_fresh(kind, t0)
                    if missing or stale:
                        self.fail(
                            f"myrun exit 0 but ROOT not fresh ({tag}); "
                            "likely stale stat from a previous myrun. "
                            f"missing={missing[:2]} stale={stale[:2]}"
                        )

            if deadline is not None and time.time() > deadline:
                snap_msg = "stat unreadable"
                if last_snap is not None:
                    snap_msg = (
                        f"done={last_snap.done_count}/"
                        f"{last_snap.expected_count} "
                        f"unfinished={last_snap.unfinished[:3]}"
                    )
                self.fail(
                    f"Recon timeout ({tag}) after {self.recon_timeout_sec}s: "
                    f"{snap_msg}"
                )

            if poll_n == 1 or poll_n % max(1, 600 // self.poll_sec) == 0:
                prog = "?"
                if last_snap is not None:
                    prog = (
                        f"{last_snap.done_count}/{last_snap.expected_count}"
                    )
                myrun_msg = "running"
                if handle is not None and self._myrun_ec_ready(handle):
                    myrun_msg = "exit 0"
                self.log(
                    f"recon ({tag}) poll #{poll_n}: progress={prog} "
                    f"myrun={myrun_msg}"
                )

            time.sleep(self.poll_sec)

        if handle is not None and not handle.foreground:
            self._wait_myrun_done(handle, block=True)

    def _wait_myrun_done(self, handle: MyrunHandle, *, block: bool) -> bool:
        if handle.foreground:
            return handle.exit_code is not None
        if handle.ec_file.is_file():
            return True
        if not block:
            return False
        wait = subprocess.run(["tmux", "wait-for", handle.token], check=False)
        if wait.returncode != 0:
            self.fail(f"tmux wait-for failed for {handle.token}")
        return True

    def run_cmd(
        self,
        cmd: list[str],
        *,
        cwd: Optional[Path] = None,
        check: bool = True,
    ) -> int:
        self.log("RUN: " + " ".join(cmd) + (f"  (cwd={cwd})" if cwd else ""))
        if self.dry_run:
            return 0
        ret = subprocess.call(cmd, cwd=str(cwd) if cwd else None)
        if check and ret != 0:
            self.fail(f"Command failed ({ret}): {' '.join(cmd)}")
        return ret

    def _tmux_available(self) -> bool:
        return bool(os.environ.get("TMUX"))

    def start_myrun(self, yml: Path, tag: str) -> Optional[MyrunHandle]:
        """Launch myrun (tmux pane or foreground). Caller polls completion."""
        yml_arg = str(yml)
        cmd_line = f"cd {self.analyzer_dir} && ./runmanager/myrun.py {yml_arg}"

        if self.dry_run:
            self.log(f"[dry-run] myrun: {cmd_line}")
            return None

        if not self._tmux_available():
            self.log(
                "[WARN] TMUX not set; running myrun in foreground "
                "(prefer: tmux new -s tpc_iter)"
            )
            ret = subprocess.call(
                ["./runmanager/myrun.py", yml_arg],
                cwd=str(self.analyzer_dir),
            )
            return MyrunHandle(
                tag=tag,
                token="",
                pane_id="",
                ec_file=Path(""),
                foreground=True,
                exit_code=ret,
            )

        token = f"tpc_recon_{tag}_{os.getpid()}_{int(time.time())}"
        ec_file = Path(f"/tmp/tpc_iter_{token}.ec")
        if ec_file.exists():
            ec_file.unlink()
        pane_script = (
            f"{cmd_line}; "
            f"ec=$?; "
            f"echo $ec > {ec_file}; "
            f"tmux wait-for -S {token}; "
            f"exit $ec"
        )
        self.log(f"tmux split-window myrun ({tag}), wait-for={token}")
        split = subprocess.run(
            [
                "tmux",
                "split-window",
                "-v",
                "-d",
                "-P",
                "-F",
                "#{pane_id}",
                pane_script,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if split.returncode != 0:
            self.fail(
                f"tmux split-window failed: {split.stderr.strip() or split.stdout}"
            )
        pane_id = split.stdout.strip()
        self.log(f"myrun pane={pane_id}")
        return MyrunHandle(
            tag=tag,
            token=token,
            pane_id=pane_id,
            ec_file=ec_file,
        )

    def cleanup_myrun(self, handle: Optional[MyrunHandle]) -> None:
        if handle is None or self.dry_run or handle.foreground:
            return
        if handle.pane_id:
            subprocess.run(
                ["tmux", "kill-pane", "-t", handle.pane_id],
                capture_output=True,
                check=False,
            )
        if handle.ec_file.is_file():
            try:
                handle.ec_file.unlink()
            except OSError:
                pass

    def recon(self, kind: str, outer_i: int, step: int) -> None:
        yml = self.yml_for(kind)
        tag = f"i{outer_i:02d}_{kind}_s{step}"
        self.notify(f"[tpc_iter] recon start {tag}\n{yml}")
        t0 = time.time()
        handle = self.start_myrun(yml, tag)
        try:
            self.wait_for_recon_complete(yml, kind, t0, handle, tag)
        finally:
            self.cleanup_myrun(handle)
        if not self.dry_run:
            snap = self._analyze_stat(
                self._read_stat_info(self.stat_json_for(yml)) or {},
                self.expected_run_keys(yml),
            )
            self.log(
                f"Recon OK ({tag}): stat {snap.done_count}/{snap.expected_count} DONE"
            )
        add = self.add_script_for(kind)
        self.run_cmd(["bash", str(add)], cwd=self.group_root)
        hadd = self.hadd_root_for(kind)
        if not self.dry_run and not hadd.is_file():
            self.fail(f"hadd output missing: {hadd}")
        self.notify(f"[tpc_iter] recon+add done {tag}\n{hadd}")

    def run_tpc(self, kind: str, cmd: str, outer_i: int, extra: list[str]) -> None:
        root = self.hadd_root_for(kind)
        run_num = self.run0_for(kind)
        mode = "trk" if kind == "trk" else "hit"
        opts = [
            sys.executable,
            str(RUN_TPC),
            str(root),
            cmd,
            "--run",
            str(run_num),
            "--mode",
            mode,
            *extra,
        ]
        self.notify(f"[tpc_iter] {cmd} start i={outer_i} kind={kind} mode={mode}")
        self.run_cmd(opts, cwd=MYANALYSIS_ROOT)
        self.notify(f"[tpc_iter] {cmd} done i={outer_i} kind={kind}")

    def offset_opts(self, kind: str) -> list[str]:
        th = "100" if kind == "hit" else "50"
        # --max-chi2-ndf 0: χ²/ndf cut off (C++ default 25 would reject many pads)
        return [
            "--threshold",
            th,
            "--max-chi2-ndf",
            "0",
            "--ave",
        ]

    def phase_opts(self, *, ofs_update: bool = False, debug: bool = False) -> list[str]:
        opts = ["--rebin", "2"]
        if ofs_update:
            opts.append("--ofs-update")
        if debug:
            opts.append("--debug")
        return opts

    def drift_opts(self, *, debug: bool = False) -> list[str]:
        return ["--debug"] if debug else []

    def run_hit_iter(self, outer_i: int) -> None:
        kind = "hit"
        step = 0

        def recon() -> None:
            nonlocal step
            step += 1
            self.recon(kind, outer_i, step=step)

        # 1) Phase + ofs-update (block head)
        recon()
        self.run_tpc(kind, "phase", outer_i, self.phase_opts(ofs_update=True))
        recon()
        self.run_tpc(kind, "offset", outer_i, self.offset_opts(kind))
        # 2) Phase (no ofs) -> offset
        recon()
        self.run_tpc(kind, "phase", outer_i, self.phase_opts())
        recon()
        self.run_tpc(kind, "offset", outer_i, self.offset_opts(kind))

    def run_trk_iter(self, outer_i: int) -> None:
        kind = "trk"
        step = 0

        def recon() -> None:
            nonlocal step
            step += 1
            self.recon(kind, outer_i, step=step)

        recon()
        self.run_tpc(kind, "offset", outer_i, self.offset_opts(kind))
        recon()
        self.run_tpc(kind, "drift", outer_i, self.drift_opts())
        # eval only (no TPCPRM write)
        self.run_tpc(kind, "phase", outer_i, self.phase_opts(debug=True))
    def run_outer_iter(self, outer_i: int) -> None:
        kind = self.kind_for_outer(outer_i)
        self.notify(
            f"[tpc_iter] OUTER START i={outer_i}/{self.max_iter} kind={kind}"
        )
        if kind == "hit":
            self.run_hit_iter(outer_i)
        else:
            self.run_trk_iter(outer_i)
        self.notify(
            f"[tpc_iter] OUTER DONE i={outer_i}/{self.max_iter} kind={kind}"
        )

    def run(self) -> None:
        if not self._tmux_available() and not self.dry_run:
            self.log(
                "[WARN] Not inside tmux. Prefer: tmux new -s tpc_iter"
            )
        self.log(
            f"Start iterate_tpc_calib max_iter={self.max_iter} "
            f"analyzer={self.analyzer_dir} dry_run={self.dry_run}"
        )
        self.log(
            f"hit runs={self.hit_runs} trk runs={self.trk_runs} "
            f"schedule=hit,hit,trk,..."
        )
        if not self.dry_run:
            self.log_dir.mkdir(parents=True, exist_ok=True)
        for i in range(1, self.max_iter + 1):
            self.run_outer_iter(i)
        self.notify(
            f"[tpc_iter] ALL DONE max_iter={self.max_iter}"
        )


def resolve_analyzer_dir() -> Path:
    # Prefer WORKDIR from hit yml if present under default analyzer path
    default = Path("/home/had/sryuta/analyzer/JPARC2025E72")
    yml = default / HIT_YML_REL
    if yml.is_file():
        with yml.open(encoding="utf-8") as f:
            data = yaml.safe_load(f)
        wd = data.get("WORKDIR")
        if wd:
            return Path(wd)
    return default


def main() -> None:
    parser = argparse.ArgumentParser(
        description="TPC calib auto iterator (hit/trk 2:1; hit: phase/offset x2; trk: offset+drift, phase --debug)"
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=6,
        help="Number of outer iterations (default 6 = two 2:1 cycles)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print plan only; do not run DST/add/calib/Discord",
    )
    parser.add_argument(
        "--analyzer-dir",
        type=Path,
        default=None,
        help="Analyzer WORKDIR (default: from hit yml WORKDIR)",
    )
    parser.add_argument(
        "--poll-sec",
        type=int,
        default=POLL_SEC,
        help="Stat JSON poll interval during recon (default 30)",
    )
    parser.add_argument(
        "--recon-timeout-sec",
        type=int,
        default=RECON_TIMEOUT_SEC,
        help="Abort recon if stat/ROOT not ready within N seconds (0=disable)",
    )
    args = parser.parse_args()
    if args.max_iter < 1:
        raise SystemExit("--max-iter must be >= 1")

    analyzer = args.analyzer_dir or resolve_analyzer_dir()
    orch = IterOrchestrator(
        analyzer_dir=analyzer,
        max_iter=args.max_iter,
        dry_run=args.dry_run,
        poll_sec=args.poll_sec,
        recon_timeout_sec=args.recon_timeout_sec,
    )
    orch.run()


if __name__ == "__main__":
    main()
