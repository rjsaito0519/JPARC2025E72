#!/usr/bin/env python3
"""
Copy selected TPCPRM fields (gain / offset) from a reference run to another run.

Usage:
  ./calibration/scripts/copy_tpcprm.py 2595 gain --ref 2594
  ./calibration/scripts/copy_tpcprm.py 2595 offset --ref 2594
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
import time
from pathlib import Path

from termcolor import colored

SCRIPT_DIR = Path(__file__).resolve().parent
MYANALYSIS_ROOT = SCRIPT_DIR.parent.parent
sys.path.append(str(MYANALYSIS_ROOT))

from lib import config
from lib.param import ParamFile

KEY_LEN = 3
FRESH_SEC = 60
ZERO_STR = "0.00000000"

TPC_PRM_DIR = config.PARAM_DIR / "TPCPRM"
E72_DIR = TPC_PRM_DIR / "e72"
EXAMPLE_PATH = TPC_PRM_DIR / "TPCParam_example"
GAIN_TEMPLATE_SRC = TPC_PRM_DIR / "TPCParam_e72_20260616"
OFFSET_TEMPLATE_SRC = E72_DIR / "TPCParam_e72_run02594"
PAD_HELPER_HH = MYANALYSIS_ROOT / "common" / "include" / "TPCPadHelper.hh"

FIELD_ATY = {"gain": 0, "offset": 2}


def run_param_path(run_num: int) -> Path:
    return E72_DIR / f"TPCParam_e72_run{run_num:05d}"


def backup_path(path: Path) -> Path:
    return path.with_name(path.name + ".bak")


def write_backup(path: Path) -> Path:
    bak = backup_path(path)
    shutil.copy2(path, bak)
    return bak


def _extract_brace_block(text: str, marker: str) -> str:
    i = text.find(marker)
    if i < 0:
        raise RuntimeError(f"marker not found in TPCPadHelper.hh: {marker}")
    j = text.find("{", i)
    if j < 0:
        raise RuntimeError(f"opening '{{' not found after {marker}")
    depth = 0
    for k in range(j, len(text)):
        if text[k] == "{":
            depth += 1
        elif text[k] == "}":
            depth -= 1
            if depth == 0:
                return text[j : k + 1]
    raise RuntimeError(f"unterminated '{{' after {marker}")


class CenterFramePads:
    """IsPadOnCenterFrame equivalent, parsed from TPCPadHelper.hh (no duplicated lists)."""

    def __init__(self, header_path: Path = PAD_HELPER_HH):
        text = header_path.read_text(encoding="utf-8", errors="replace")
        pad_block = _extract_brace_block(text, "padParameter[NumOfLayersTPC][NPadParameter]")
        rows = re.findall(r"\{([^{}]+)\}", pad_block)
        if len(rows) != 32:
            raise RuntimeError(f"expected 32 padParameter rows, got {len(rows)}")
        self.n_pad: list[int] = []
        for row in rows:
            parts = [p.strip() for p in row.split(",") if p.strip()]
            if len(parts) < 2:
                raise RuntimeError(f"bad padParameter row: {row}")
            self.n_pad.append(int(float(parts[1])))

        cf_block = _extract_brace_block(text, "padOnCenterFrame[]")
        self.pad_ids = {int(x) for x in re.findall(r"\d+", cf_block)}
        if not self.pad_ids:
            raise RuntimeError("padOnCenterFrame is empty")

    def get_pad_id(self, layer: int, row: int) -> int:
        if row < 0:
            return -1
        pad_id = 0
        for layer_i in range(layer):
            pad_id += self.n_pad[layer_i]
        return pad_id + row

    def is_on_center_frame(self, layer: int, row: int) -> bool:
        return self.get_pad_id(layer, row) in self.pad_ids


def filter_by_aty(data: dict[str, list[str]], aty: int) -> dict[str, list[str]]:
    aty_s = str(aty)
    out: dict[str, list[str]] = {}
    for key, vals in data.items():
        parts = key.split("-")
        if len(parts) == 3 and parts[2] == aty_s:
            out[key] = vals
    return out


def overlay_aty(dest: ParamFile, src: ParamFile, aty: int) -> tuple[int, int, int]:
    """Copy src ATY rows onto matching dest keys. Returns (updated, src_not_in_dest, dest_not_in_src)."""
    src_aty = filter_by_aty(src.get_data_dict(key_len=KEY_LEN), aty)
    dest_all = dest.get_data_dict(key_len=KEY_LEN)
    dest_aty = filter_by_aty(dest_all, aty)
    missing_in_dest = sum(1 for k in src_aty if k not in dest_all)
    missing_in_src = sum(1 for k in dest_aty if k not in src_aty)
    n_updated = dest.update(src_aty, key_len=KEY_LEN, start_col=KEY_LEN)
    return n_updated, missing_in_dest, missing_in_src


def zero_center_frame_gain(param: ParamFile, pads: CenterFramePads) -> tuple[int, int, int]:
    """Force ATY=0 p0/p1 = 0 on center-frame pads. Returns (set, already_zero, n_on_frame)."""
    n_set = 0
    n_already = 0
    n_on_frame = 0
    for i, (is_comment, parts) in enumerate(param.records):
        if is_comment or len(parts) < 5:
            continue
        try:
            layer = int(parts[0])
            row = int(parts[1])
            aty = int(parts[2])
        except ValueError:
            continue
        if aty != 0:
            continue
        if not pads.is_on_center_frame(layer, row):
            continue
        n_on_frame += 1
        try:
            p0 = float(parts[3])
            p1 = float(parts[4])
        except ValueError:
            p0, p1 = 1.0, 1.0
        if p0 == 0.0 and p1 == 0.0:
            n_already += 1
            continue
        new_parts = list(parts)
        new_parts[3] = ZERO_STR
        new_parts[4] = ZERO_STR
        param.records[i] = (False, new_parts)
        n_set += 1
    return n_set, n_already, n_on_frame


def confirm_yes() -> bool:
    if not sys.stdin.isatty() or not sys.stdout.isatty():
        print(
            colored(
                "[Error] 対話確認が必要なので、ターミナル (TTY) で実行してください。",
                "red",
            )
        )
        sys.exit(1)
    ans = input("実行しますか? [y/N] ").strip().lower()
    return ans in ("y", "yes")


def warn_if_fresh(path: Path) -> None:
    age = time.time() - path.stat().st_mtime
    if age < FRESH_SEC:
        print(
            colored(
                f"[WARN] {path.name} は {age:.1f} 秒前に更新されています。"
                f"{FRESH_SEC} 秒以内の上書きは、直前の較正を潰している可能性があります。",
                "yellow",
            )
        )


def init_example_template(pads: CenterFramePads | None = None) -> None:
    """One-shot: overlay gain(20260616) + offset(run02594) onto TPCParam_example, zero CF gain."""
    if not EXAMPLE_PATH.is_file():
        raise SystemExit(f"[Error] template not found: {EXAMPLE_PATH}")
    if not GAIN_TEMPLATE_SRC.is_file():
        raise SystemExit(f"[Error] gain source not found: {GAIN_TEMPLATE_SRC}")
    if not OFFSET_TEMPLATE_SRC.is_file():
        raise SystemExit(f"[Error] offset source not found: {OFFSET_TEMPLATE_SRC}")
    if pads is None:
        pads = CenterFramePads()

    bak = write_backup(EXAMPLE_PATH)
    dest = ParamFile(EXAMPLE_PATH)
    gain_src = ParamFile(GAIN_TEMPLATE_SRC)
    ofs_src = ParamFile(OFFSET_TEMPLATE_SRC)

    n_gain, miss_g_dest, miss_g_src = overlay_aty(dest, gain_src, FIELD_ATY["gain"])
    n_ofs, miss_o_dest, miss_o_src = overlay_aty(dest, ofs_src, FIELD_ATY["offset"])
    n_set, n_already, n_cf = zero_center_frame_gain(dest, pads)
    dest.write()

    print(colored(f"[INFO] backup: {bak}", "cyan"))
    print(colored(f"[INFO] gain overlay (ATY=0): {n_gain} rows from {GAIN_TEMPLATE_SRC.name}", "green"))
    if miss_g_dest or miss_g_src:
        print(colored(f"[WARN] gain key mismatch: src-only={miss_g_dest} dest-only={miss_g_src}", "yellow"))
    print(colored(f"[INFO] offset overlay (ATY=2): {n_ofs} rows from {OFFSET_TEMPLATE_SRC.name}", "green"))
    if miss_o_dest or miss_o_src:
        print(colored(f"[WARN] offset key mismatch: src-only={miss_o_dest} dest-only={miss_o_src}", "yellow"))
    print(
        colored(
            f"[INFO] center-frame ATY=0: on_frame={n_cf} already_zero={n_already} set_to_zero={n_set}",
            "green",
        )
    )
    print(colored(f"[INFO] wrote {EXAMPLE_PATH}", "green"))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Copy TPCPRM gain (ATY=0) or offset (ATY=2, p0+p1) from --ref run to dest run."
    )
    parser.add_argument("run_num", type=int, help="Destination run number")
    parser.add_argument(
        "field",
        choices=["gain", "offset"],
        help="gain: ATY=0; offset: ATY=2 (time offset p0 and vdrift p1)",
    )
    parser.add_argument(
        "--ref",
        type=int,
        required=True,
        metavar="RUN",
        help="Source run number to copy from",
    )
    args = parser.parse_args()

    src_path = run_param_path(args.ref)
    dest_path = run_param_path(args.run_num)
    aty = FIELD_ATY[args.field]

    if not src_path.is_file():
        print(colored(f"[Error] コピー元がありません: {src_path}", "red"))
        sys.exit(1)
    if not EXAMPLE_PATH.is_file():
        print(colored(f"[Error] テンプレがありません: {EXAMPLE_PATH}", "red"))
        sys.exit(1)

    dest_exists = dest_path.is_file()
    pads = CenterFramePads() if args.field == "gain" else None

    print(colored(">>> TPCPRM copy plan <<<", "cyan", attrs=["bold"]))
    print(f"  field : {args.field} (ATY={aty})")
    print(f"  src   : {src_path}")
    print(f"  dest  : {dest_path}")
    if dest_exists:
        print(f"  dest  : 既存ファイルを {dest_path.name}.bak に退避してから上書き")
        warn_if_fresh(dest_path)
    else:
        print(f"  dest  : 無いので {EXAMPLE_PATH.name} をコピーしてから該当 ATY だけ上書き")

    src = ParamFile(src_path)
    src_n = len(filter_by_aty(src.get_data_dict(key_len=KEY_LEN), aty))
    print(f"  src ATY={aty} rows: {src_n}")

    if not confirm_yes():
        print(colored("[SKIP] 中止しました。", "yellow"))
        sys.exit(0)

    if dest_exists:
        bak = write_backup(dest_path)
        print(colored(f"[INFO] backup: {bak}", "cyan"))
        dest = ParamFile(dest_path)
    else:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(EXAMPLE_PATH, dest_path)
        print(colored(f"[INFO] created {dest_path.name} from {EXAMPLE_PATH.name}", "green"))
        dest = ParamFile(dest_path)

    n_upd, miss_dest, miss_src = overlay_aty(dest, src, aty)
    if miss_dest or miss_src:
        print(
            colored(
                f"[WARN] key mismatch ATY={aty}: src-only={miss_dest} dest-only={miss_src}",
                "yellow",
            )
        )
    print(colored(f"[INFO] overlay ATY={aty}: {n_upd} rows", "green"))

    if args.field == "gain":
        assert pads is not None
        n_set, n_already, n_cf = zero_center_frame_gain(dest, pads)
        print(
            colored(
                f"[INFO] center-frame ATY=0: on_frame={n_cf} already_zero={n_already} set_to_zero={n_set}",
                "green",
            )
        )

    dest.write()
    print(colored(f"[SUCCESS] wrote {dest_path}", "green", attrs=["bold"]))


if __name__ == "__main__":
    main()
