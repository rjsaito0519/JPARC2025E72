#!/usr/bin/env python3
"""Create DST runlist / conf / symlink for TPCHitBcOut, TPCHelix, TPCTracking."""

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config

SCRATCH_RUNLIST_DIR = "scratch"
SUB_DIR = config.SUB_DIR

# Standard legacy flat names -> new TAG
MIGRATE_PATTERNS = [
    (re.compile(r"^dst_tpchit_bcout_tracking_run(\d{5})\.root$"), "TPCHitBcOut"),
    (re.compile(r"^dst_tpchelix_run(\d{5})\.root$"), "TPCHelix"),
    (re.compile(r"^dst_tpctracking_run(\d{5})\.root$"), "TPCTracking"),
    # mis-tagged from early migrate
    (re.compile(r"^run(\d{5})_TPCTrack\.root$"), "TPCTracking"),
]

MODES = {
    "tpchit_bcout": {
        "cli": "--tpchit-bcout",
        "flag": "tpchit_bcout",
        "bin": "./bin/DstTPCHitBcOutTracking",
        "tag": "TPCHitBcOut",
        "mode_label": "tpchit_bcout",
        "yml_prefix": "dst_tpchit_bcout",
        "tpl_yml": "dst_tpchit_bcout_example.yml",
        "tpl_conf": "analyzer_e72_dst_tpchit_bcout_example.conf",
        "dstin": ["TPC", "BcOut"],
        "unit": 50000,
        "option_key": "option",
        "option_val": "-n 2",
        "required_inputs": ["TPC", "BcOut"],
    },
    "tpchelix": {
        "cli": "--tpchelix",
        "flag": "tpchelix",
        "bin": "./bin/DstTPCHelixTracking",
        "tag": "TPCHelix",
        "mode_label": "tpchelix",
        "yml_prefix": "dst_tpchelix",
        "tpl_yml": "dst_tpchelix_example.yml",
        "tpl_conf": "analyzer_e72_dst_tpchelix_example.conf",
        "dstin": ["TPC"],
        "unit": 3000,
        "option_key": "moption",
        "option_val": "-n 4",
        "required_inputs": ["TPC"],
    },
    "tpctrack": {
        "cli": "--tpctrack",
        "flag": "tpctrack",
        "bin": "./bin/DstTPCTracking",
        "tag": "TPCTracking",
        "mode_label": "tpctrack",
        "yml_prefix": "dst_tpctrack",
        "tpl_yml": "dst_tpctrack_example.yml",
        "tpl_conf": "analyzer_e72_dst_tpctrack_example.conf",
        "dstin": ["TPC"],
        "unit": 25000,
        "option_key": "option",
        "option_val": "-n 2",
        "required_inputs": ["TPC"],
    },
}

DEFAULTS = {
    "USER:": "param/USER/UserParam_e72_20260403",
    "TPCPRM:": "param/TPCPRM/TPCParam_example",
    "TPCPOS:": "param/TPCPOS/TPCPositionCorrectionMap_KinFit_th50_1",
    "DCGEO:": "param/DCGEO/DCGeomParam_e72_example",
    "FLDMAP:": "param/FLDMAP/ShsFieldMap_20210526_Extrapolated",
}

HSFLD_DEFAULTS = {
    "HSFLDCALIB:": "0.94838",
    "HSFLDCALC:": "0.90607",
    "HSFLDHALL:": "0.90",
}

FILE_KEYS = ("USER:", "TPCPRM:", "TPCPOS:", "DCGEO:", "FLDMAP:")


def colored(text, color):
    codes = {"green": "32", "cyan": "36", "yellow": "33", "red": "31"}
    return f"\033[{codes.get(color, '0')}m{text}\033[0m"


def param_abs_path(rel_path: str) -> Path:
    if rel_path.startswith("param/"):
        return config.PARAM_DIR / rel_path[len("param/"):]
    return config.ANALYZER_DIR / rel_path


def ensure_data_symlink(link_path: Path, target: Path):
    try:
        if link_path.is_symlink() or link_path.exists() or link_path.is_file():
            link_path.unlink()
        os.symlink(target.resolve(), link_path)
    except OSError as e:
        print(colored(f"Warning: Failed to create symlink: {e}", "yellow"))


def input_exists(path: Path) -> bool:
    if path.is_symlink():
        try:
            return path.resolve().is_file()
        except OSError:
            return False
    return path.is_file()


def check_required_inputs(run_num: int, mode: dict) -> list:
    missing = []
    run_tag = f"run{run_num:05d}"
    for key in mode["required_inputs"]:
        path = config.DATA_DIR / f"{run_tag}_{key}.root"
        if not input_exists(path):
            missing.append(path)
    return missing


def read_conf_map(conf_path: Path) -> dict:
    result = {}
    if not conf_path.exists():
        return result
    for line in conf_path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and not line.lstrip().startswith("#"):
            result[parts[0]] = parts[1]
    return result


def read_digit_cmap(conf_path: Path):
    digit = cmap = None
    for line in conf_path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            if parts[0] == "DIGIT:":
                digit = parts[1]
            elif parts[0] == "CMAP:":
                cmap = parts[1]
    if digit is None or cmap is None:
        print(colored(f"Error: DIGIT/CMAP not found in {conf_path}", "red"))
        sys.exit(1)
    return digit, cmap


def era_conf_name(run_num: int) -> str:
    if run_num < 2568:
        return "analyzer_e72_example1.conf"
    if run_num < 3642:
        return "analyzer_e72_example2.conf"
    return "analyzer_e72_example3.conf"


def find_run_specific(key: str, run_num: int):
    """Return relative param/ path if a run-specific file exists."""
    n = f"{run_num:05d}"
    candidates = []
    if key == "TPCPRM:":
        candidates = [
            f"param/TPCPRM/{SUB_DIR}/TPCParam_e72_run{n}",
            f"param/TPCPRM/TPCParam_e72_run{n}",
            f"param/TPCPRM/{SUB_DIR}/TPCParam_run{n}",
        ]
    elif key == "USER:":
        candidates = [
            f"param/USER/{SUB_DIR}/UserParam_run{n}",
            f"param/USER/UserParam_run{n}",
        ]
    elif key == "DCGEO:":
        candidates = [
            f"param/DCGEO/{SUB_DIR}/DCGeomParam_run{n}_K",
            f"param/DCGEO/{SUB_DIR}/DCGeomParam_run{n}_Pi",
            f"param/DCGEO/{SUB_DIR}/DCGeomParam_run{n}",
            f"param/DCGEO/DCGeomParam_run{n}",
        ]
    elif key == "TPCPOS:":
        candidates = [
            f"param/TPCPOS/{SUB_DIR}/TPCPositionCorrectionMap_run{n}",
        ]
    for rel in candidates:
        if param_abs_path(rel).exists():
            return rel
    return None


def find_ref_conf(param_run: int, mode_label: str):
    conf_dir = config.PARAM_DIR / "conf" / SUB_DIR
    primary = conf_dir / f"analyzer_run{param_run:05d}_dst_{mode_label}.conf"
    if primary.exists():
        return primary, False
    for pat in [
        f"analyzer_run{param_run:05d}_dst_{mode_label}_*.conf",
        f"analyzer_run{param_run:05d}_bcout_*.conf",
    ]:
        cands = sorted(conf_dir.glob(pat))
        if cands:
            return cands[0], True
    for name in (
        f"analyzer_run{param_run:05d}_tpchit.conf",
        f"analyzer_run{param_run:05d}_bcout_K.conf",
        f"analyzer_run{param_run:05d}_bcout_Pi.conf",
    ):
        p = conf_dir / name
        if p.exists():
            return p, True
    return None, False


def resolve_key_path(key: str, run_num: int, mode: dict, ref):
    # Tier 1: run-specific file for this run
    own = find_run_specific(key, run_num)
    if own is not None:
        return own

    # Tier 2: --ref
    if ref is not None:
        ref_own = find_run_specific(key, ref)
        if ref_own is not None:
            print(colored(
                f"Warning: --ref {ref} using run-specific {key} {ref_own}",
                "yellow",
            ))
            return ref_own
        ref_conf, is_fb = find_ref_conf(ref, mode["mode_label"])
        if ref_conf is not None:
            cmap = read_conf_map(ref_conf)
            if key in cmap and param_abs_path(cmap[key]).exists():
                note = ref_conf.name + (" (fallback conf)" if is_fb else "")
                print(colored(
                    f"Warning: --ref {ref} {key} from {note}: {cmap[key]}",
                    "yellow",
                ))
                return cmap[key]

    # Tier 3: template default
    return DEFAULTS[key]


def validate_dcgeo(dcgeo_rel: str) -> None:
    path = param_abs_path(dcgeo_rel)
    if not path.exists():
        print(colored(f"Error: DCGEO not found: {path}", "red"))
        sys.exit(1)
    text = path.read_text()
    has_tpc = bool(re.search(r"(?m)^\s*4001\b", text))
    has_tgt = bool(re.search(r"(?m)^\s*4002\b", text))
    if not has_tpc or not has_tgt:
        print(colored(
            f"Error: DCGEO missing HypTPC(4001)/Target(4002): {path}",
            "red",
        ))
        sys.exit(1)


def validate_conf(conf_path: Path) -> None:
    cmap = read_conf_map(conf_path)
    for key in ("HSFLDCALIB:", "HSFLDCALC:", "HSFLDHALL:"):
        if key not in cmap:
            print(colored(f"Error: {key} missing in {conf_path}", "red"))
            sys.exit(1)
        try:
            float(cmap[key])
        except ValueError:
            print(colored(f"Error: {key} not numeric in {conf_path}: {cmap[key]}", "red"))
            sys.exit(1)

    for key in FILE_KEYS:
        if key not in cmap:
            print(colored(f"Error: {key} missing in {conf_path}", "red"))
            sys.exit(1)
        abs_p = param_abs_path(cmap[key])
        if not abs_p.exists():
            print(colored(f"Error: {key} file not found: {abs_p}", "red"))
            sys.exit(1)

    validate_dcgeo(cmap["DCGEO:"])


def ensure_dst_conf(run_num: int, mode: dict, ref) -> str:
    conf_dir = config.PARAM_DIR / "conf"
    sub_conf = conf_dir / SUB_DIR
    sub_conf.mkdir(parents=True, exist_ok=True)

    tpl_path = conf_dir / mode["tpl_conf"]
    era_path = conf_dir / era_conf_name(run_num)
    digit, cmap_era = read_digit_cmap(era_path)

    resolved = {k: resolve_key_path(k, run_num, mode, ref) for k in DEFAULTS}
    # Prefer HSFLD from --ref conf if present
    hsfld = dict(HSFLD_DEFAULTS)
    if ref is not None:
        ref_conf, _ = find_ref_conf(ref, mode["mode_label"])
        if ref_conf is not None:
            ref_map = read_conf_map(ref_conf)
            for k in HSFLD_DEFAULTS:
                if k in ref_map:
                    hsfld[k] = ref_map[k]

    out_lines = []
    seen_hsfld = set()
    for line in tpl_path.read_text().splitlines(keepends=True):
        stripped = line.lstrip()
        if stripped.startswith("#") or not stripped.strip():
            out_lines.append(line)
            continue
        parts = stripped.split()
        if not parts:
            out_lines.append(line)
            continue
        key = parts[0]
        prefix = line[: len(line) - len(stripped)]
        newline = "\n" if line.endswith("\n") else ""
        if key == "DIGIT:":
            out_lines.append(f"{prefix}DIGIT:\t\t{digit}{newline}")
        elif key == "CMAP:":
            out_lines.append(f"{prefix}CMAP:\t\t{cmap_era}{newline}")
        elif key in resolved:
            out_lines.append(f"{prefix}{key}\t\t{resolved[key]}{newline}")
        elif key in hsfld:
            out_lines.append(f"{prefix}{key}\t{hsfld[key]}{newline}")
            seen_hsfld.add(key)
        else:
            out_lines.append(line)

    for key, val in hsfld.items():
        if key not in seen_hsfld:
            out_lines.append(f"{key}\t{val}\n")

    target = sub_conf / f"analyzer_run{run_num:05d}_dst_{mode['mode_label']}.conf"
    target.write_text("".join(out_lines))
    validate_conf(target)
    return f"param/conf/{SUB_DIR}/{target.name}"


def write_dst_runlist(path: Path, mode: dict, runs_info: list):
    runlist_dir = config.ANALYZER_DIR / "runmanager/runlist"
    tpl_path = runlist_dir / mode["tpl_yml"]
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(tpl_path) as f_tpl:
        tpl_lines = f_tpl.readlines()

    with open(path, "w") as f_out:
        for line in tpl_lines:
            if line.strip() == "RUN:":
                break
            f_out.write(line)
        f_out.write("RUN:\n")
        for info in runs_info:
            f_out.write(f"  {info['label']}\n")
            f_out.write(f"    bin: {info['bin']}\n")
            f_out.write(f"    conf: {info['conf']}\n")
            f_out.write(f"    data: {info['data']}\n")
            f_out.write(f"    root: {info['root']}\n")
            f_out.write(f"    dstin_root: {info['dstin_root']}\n")
            f_out.write("    dstin:\n")
            for d in info["dstin"]:
                f_out.write(f"      - {d}\n")
            f_out.write(f"    unit: {info['unit']}\n")
            f_out.write(f"    {info['option_key']}: \"{info['option_val']}\"\n")


def batch_runlist_name(mode: dict, run_nums: list) -> str:
    run_nums = sorted(run_nums)
    prefix = mode["yml_prefix"]
    if len(run_nums) == 1:
        return f"{prefix}_run{run_nums[0]:05d}.yml"
    return f"{prefix}_{run_nums[0]:05d}-{run_nums[-1]:05d}_{len(run_nums)}runs.yml"


def setup_dst_run(run_num: int, mode: dict, ref):
    run_tag = f"run{run_num:05d}"
    tag = mode["tag"]

    missing = check_required_inputs(run_num, mode)
    if missing:
        for p in missing:
            print(colored(
                f"Error: missing input for run {run_num:05d} {mode['cli']}: {p}\n"
                f"       (expected symlink or file; create via create_runlist.py "
                f"--tpchit / --bcout first)",
                "red",
            ))
        return None

    (config.DECODE_DIR / run_tag).mkdir(parents=True, exist_ok=True)
    (config.SCRATCH_DIR / run_tag).mkdir(parents=True, exist_ok=True)

    root_name = f"{run_tag}_{tag}.root"
    actual_root = config.DECODE_DIR / run_tag / root_name
    symlink_path = config.DATA_DIR / root_name
    ensure_data_symlink(symlink_path, actual_root)

    conf_rel = ensure_dst_conf(run_num, mode, ref)

    run_info = {
        "run_num": run_num,
        "label": f"{run_tag}_{mode['mode_label']}:",
        "bin": mode["bin"],
        "conf": conf_rel,
        "data": f"rawdata/{run_tag}.dat",
        "root": str(actual_root),
        "dstin_root": str(config.DATA_DIR),
        "dstin": list(mode["dstin"]),
        "unit": mode["unit"],
        "option_key": mode["option_key"],
        "option_val": mode["option_val"],
    }

    runlist_path = (
        config.ANALYZER_DIR / "runmanager/runlist" / SUB_DIR
        / f"{mode['yml_prefix']}_{run_tag}.yml"
    )
    write_dst_runlist(runlist_path, mode, [run_info])

    return runlist_path, symlink_path, actual_root, run_info


def classify_dst_entry(name: str):
    for regex, tag in MIGRATE_PATTERNS:
        m = regex.match(name)
        if m:
            return int(m.group(1)), tag, "standard"
    if name.startswith("dst_") or name.startswith("bk_dst_"):
        return None, None, "irregular"
    return None, None, "ignore"


def migrate_one(run_num: int, tag: str, src: Path, dry_run: bool) -> bool:
    run_tag = f"run{run_num:05d}"
    dest_dir = config.DECODE_DIR / run_tag
    dest_file = dest_dir / f"{run_tag}_{tag}.root"
    symlink_path = config.DATA_DIR / f"{run_tag}_{tag}.root"

    if dest_file.exists() and dest_file.is_file():
        print(colored(f"[OK] run {run_num:05d} {tag}: already migrated", "green"))
        if not dry_run:
            ensure_data_symlink(symlink_path, dest_file)
            if src.exists() or src.is_symlink():
                if src.resolve() != dest_file.resolve():
                    # leave old only if different; if same inode after move skip
                    pass
        return True

    print(colored(f"[MOVE] {src.name} -> {dest_file}", "cyan"))
    if dry_run:
        return True

    dest_dir.mkdir(parents=True, exist_ok=True)
    if dest_file.exists() or dest_file.is_symlink():
        dest_file.unlink()
    shutil.move(str(src), str(dest_file))
    ensure_data_symlink(symlink_path, dest_file)
    return True


def migrate_dst_roots(run_nums=None, dry_run: bool = False):
    data_dir = config.DATA_DIR
    irregular = []
    to_move = []  # (run_num, tag, path)

    for entry in sorted(data_dir.iterdir()):
        if not entry.name.endswith(".root"):
            continue
        run_num, tag, kind = classify_dst_entry(entry.name)
        if kind == "irregular":
            irregular.append(entry.name)
            continue
        if kind != "standard":
            continue
        if run_nums is not None and run_num not in run_nums:
            continue
        to_move.append((run_num, tag, entry))

    if irregular:
        print(colored("\n--- Irregular DST filenames (left untouched) ---", "yellow"))
        for name in irregular:
            print(f"  {name}")

    # detect collisions: multiple sources -> same (run, tag)
    by_key = {}
    for run_num, tag, path in to_move:
        key = (run_num, tag)
        by_key.setdefault(key, []).append(path)

    migrated = 0
    for (run_num, tag), paths in sorted(by_key.items()):
        if len(paths) > 1:
            print(colored(
                f"[SKIP] run {run_num:05d} {tag}: multiple sources — consult manually",
                "yellow",
            ))
            for p in paths:
                print(f"         {p}")
            continue
        if migrate_one(run_num, tag, paths[0], dry_run=dry_run):
            migrated += 1

    print(colored(f"\n[MIGRATE] processed {migrated} standard DST file(s)", "green"))
    return migrated


def main():
    parser = argparse.ArgumentParser(
        prog="dst_create_runlist",
        description="Create DST conf/runlist/symlink for TPCHitBcOut / TPCHelix / TPCTracking.",
    )
    parser.add_argument("run_nums", type=int, nargs="*", help="Run number(s)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--tpchit-bcout", action="store_true", help="DstTPCHitBcOutTracking")
    group.add_argument("--tpchelix", action="store_true", help="DstTPCHelixTracking")
    group.add_argument("--tpctrack", action="store_true", help="DstTPCTracking (tag: TPCTracking)")
    parser.add_argument("--ref", type=int, default=None, help="Borrow params from this run")
    parser.add_argument(
        "--migrate", action="store_true",
        help="Migrate standard flat dst_*.root to root/runNNNNN/runNNNNN_{TAG}.root",
    )
    parser.add_argument("--dry-run", action="store_true", help="Show migrate actions only")
    args = parser.parse_args()

    mode = None
    if args.tpchit_bcout:
        mode = MODES["tpchit_bcout"]
    elif args.tpchelix:
        mode = MODES["tpchelix"]
    elif args.tpctrack:
        mode = MODES["tpctrack"]

    if args.migrate:
        migrate_dst_roots(args.run_nums or None, dry_run=args.dry_run)
        if args.dry_run or mode is None:
            return

    if mode is None:
        if not args.migrate:
            print(colored(
                "Error: choose one of --tpchit-bcout / --tpchelix / --tpctrack "
                "(or --migrate only)",
                "red",
            ))
            sys.exit(1)
        return

    if not args.run_nums:
        print(colored("Error: at least one run number is required", "red"))
        sys.exit(1)

    run_nums = sorted(set(args.run_nums))
    all_infos = []
    failed = []

    for run_num in run_nums:
        result = setup_dst_run(run_num, mode, args.ref)
        if result is None:
            failed.append(run_num)
            continue
        yml, link, root, run_info = result
        all_infos.append(run_info)
        print(colored(f"\n[SUCCESS] DST {mode['tag']} setup for run {run_num:05d}", "green"))
        print(f"  - Conf: {run_info['conf']}")
        print(f"  - Root: {root}")
        print(f"  - Link: {link}")
        print(colored(f"  - Runlist: {yml.relative_to(config.ANALYZER_DIR)}", "cyan"))
        print(colored("\n--- Runlist File Content ---", "cyan"))
        with open(yml) as f:
            print(f.read())

    if len(all_infos) > 1:
        batch_name = batch_runlist_name(mode, [x["run_num"] for x in all_infos])
        batch_yml = (
            config.ANALYZER_DIR / "runmanager/runlist" / SCRATCH_RUNLIST_DIR / batch_name
        )
        write_dst_runlist(batch_yml, mode, all_infos)
        print(colored(
            f"\n[BATCH] Combined runlist ({len(all_infos)} runs)",
            "green",
        ))
        print(colored(f"  - {batch_yml.relative_to(config.ANALYZER_DIR)}", "cyan"))

    if failed:
        print(colored(
            f"\nError: failed for run(s): {', '.join(f'{r:05d}' for r in failed)}",
            "red",
        ))
        sys.exit(1)


if __name__ == "__main__":
    main()
