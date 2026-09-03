#!/usr/bin/env python3

# ---------------------------------------------------------------------------
import argparse
import re
import sys
import os
import shutil
from pathlib import Path

# Add shared library path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config

TPCHIT_CONF_TEMPLATE = "analyzer_e72_tpchit_example.conf"
TPCHIT_USER = config.PARAM_DIR / "USER" / "UserParam_e72_tpc_0"
TPCHIT_HDPRM = config.PARAM_DIR / "HDPRM" / "HodoParam_02293"
TPCHIT_SCRATCH_RUNLIST_DIR = "scratch"

TpcRunRe = re.compile(r"^tpc_run(\d{5})\.root$")
RunTpcRe = re.compile(r"^run(\d{5})_TPC\.root$")
IrregularTpcRe = re.compile(r"^tpc_run(\d{5})_.+\.root$")


def colored(text, color):
    codes = {"green": "32", "cyan": "36", "yellow": "33", "red": "31"}
    return f"\033[{codes.get(color, '0')}m{text}\033[0m"


def set_central_momentum(path, mom):
    lines = path.read_text().splitlines(keepends=True)
    found = False
    out = []
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("CentralMomentum:") and not stripped.startswith("#"):
            prefix_ws = line[: len(line) - len(stripped)]
            newline = "\n" if line.endswith("\n") else ""
            out.append(f"{prefix_ws}CentralMomentum: {mom}{newline}")
            found = True
        else:
            out.append(line)
    if not found:
        print(colored(f"Error: CentralMomentum: not found in {path}", "red"))
        sys.exit(1)
    path.write_text("".join(out))


def check_cobo_consistency():
    if not TPCHIT_USER.exists() or not TPCHIT_HDPRM.exists():
        print(colored("Warning: TPCHit USER/HDPRM not found; skip CoBo check", "yellow"))
        return

    cobo_min = cobo_max = None
    for line in TPCHIT_USER.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0] == "COBO_TDC":
            cobo_min = float(parts[1])
            cobo_max = float(parts[2])
            break
    if cobo_min is None:
        print(colored("Warning: COBO_TDC not found in UserParam_e72_tpc_0", "yellow"))
        return

    hdprm_offsets = []
    in_cobo = False
    for line in TPCHIT_HDPRM.read_text().splitlines():
        if line.startswith("#COBO"):
            in_cobo = True
            continue
        if in_cobo:
            if line.startswith("#") or not line.strip():
                break
            parts = line.split()
            if len(parts) >= 6:
                hdprm_offsets.append(float(parts[5]))

    if not hdprm_offsets:
        print(colored("Warning: #COBO section not found in HodoParam_02293", "yellow"))
        return

    offset = hdprm_offsets[0]
    if offset < cobo_min or offset > cobo_max:
        print(colored(
            f"Warning: CoBo mismatch — HDPRM offset {offset:.0f} "
            f"outside USER COBO_TDC [{cobo_min:.0f}, {cobo_max:.0f}]",
            "yellow",
        ))
    else:
        print(colored(
            f"CoBo OK: HDPRM offset {offset:.0f} within "
            f"USER COBO_TDC [{cobo_min:.0f}, {cobo_max:.0f}]",
            "green",
        ))


def resolve_actual_file(path: Path) -> Path:
    if path.is_symlink():
        target = path.resolve()
        if target.exists():
            return target
    if path.is_file() and not path.is_symlink():
        return path.resolve()
    return path.resolve()


def ensure_data_symlink(link_path: Path, target: Path):
    """Replace existing symlink or file at link_path with a new symlink."""
    try:
        if link_path.is_symlink() or link_path.exists() or link_path.is_file():
            link_path.unlink()
        os.symlink(target.resolve(), link_path)
    except OSError as e:
        print(colored(f"Warning: Failed to create symlink: {e}", "yellow"))


def migrate_tpchit_root(run_num: int, dry_run: bool = False):
    data_dir = config.DATA_DIR
    decode_dir = config.DECODE_DIR
    run_tag = f"run{run_num:05d}"
    dest_dir = decode_dir / run_tag
    dest_file = dest_dir / f"{run_tag}_TPC.root"
    symlink_path = data_dir / f"{run_tag}_TPC.root"
    legacy_path = data_dir / f"tpc_{run_tag}.root"

    sources = []
    if legacy_path.exists() or legacy_path.is_symlink():
        sources.append(resolve_actual_file(legacy_path))
    if symlink_path.exists() or symlink_path.is_symlink():
        resolved = resolve_actual_file(symlink_path)
        if resolved not in sources:
            sources.append(resolved)

    real_files = [p for p in sources if p.is_file()]
    if len(real_files) > 1:
        print(colored(
            f"[SKIP] run {run_num:05d}: multiple real files — consult manually",
            "yellow",
        ))
        for p in real_files:
            print(f"         {p}")
        return False

    if dest_file.exists() and dest_file.is_file():
        if symlink_path.is_symlink() and symlink_path.resolve() == dest_file.resolve():
            print(colored(f"[OK] run {run_num:05d}: already migrated", "green"))
            return True
        if not real_files:
            print(colored(f"[OK] run {run_num:05d}: dest exists, fixing symlink", "green"))
            if not dry_run:
                dest_dir.mkdir(parents=True, exist_ok=True)
                if symlink_path.exists() or symlink_path.is_symlink():
                    symlink_path.unlink()
                os.symlink(dest_file.resolve(), symlink_path)
                if legacy_path.exists() or legacy_path.is_symlink():
                    legacy_path.unlink()
            return True

    if not real_files:
        print(colored(f"[SKIP] run {run_num:05d}: no source ROOT found", "yellow"))
        return False

    src = real_files[0]
    print(colored(f"[MOVE] run {run_num:05d}: {src.name} -> {dest_file}", "cyan"))
    if dry_run:
        return True

    dest_dir.mkdir(parents=True, exist_ok=True)
    if dest_file.exists() or dest_file.is_symlink():
        dest_file.unlink()
    shutil.move(str(src), str(dest_file))

    for p in (symlink_path, legacy_path):
        if p.exists() or p.is_symlink():
            p.unlink()
    os.symlink(dest_file.resolve(), symlink_path)
    return True


def migrate_all_tpchit_roots(run_nums=None, dry_run: bool=False):
    data_dir = config.DATA_DIR
    consult = []

    for entry in sorted(data_dir.iterdir()):
        if not entry.name.endswith(".root"):
            continue
        m_irreg = IrregularTpcRe.match(entry.name)
        if m_irreg:
            consult.append(entry.name)
            continue

    if consult:
        print(colored("\n--- Irregular filenames (manual review) ---", "yellow"))
        for name in consult:
            print(f"  {name}")

    if run_nums is None:
        run_set = set()
        for entry in data_dir.iterdir():
            m1 = TpcRunRe.match(entry.name)
            m2 = RunTpcRe.match(entry.name)
            if m1:
                run_set.add(int(m1.group(1)))
            elif m2:
                run_set.add(int(m2.group(1)))
        run_nums = sorted(run_set)

    migrated = 0
    for run_num in run_nums:
        if migrate_tpchit_root(run_num, dry_run=dry_run):
            migrated += 1

    print(colored(f"\n[MIGRATE] processed {migrated}/{len(run_nums)} run(s)", "green"))
    return migrated


RunLinkRe = re.compile(r"^(run\d{5})_(Hodo|BcOut|BcIn|D5)\.root$")
MiscRunRe = re.compile(r"^(run\d{5})_.+\.root$")


def find_calibration_final(run_tag: str, kind: str):
    decode_dir = config.DECODE_DIR / run_tag
    if not decode_dir.exists():
        return None
    patterns = {
        "Hodo": f"{run_tag}_Hodo_hdphc_*.root",
        "BcOut": f"{run_tag}_BcOut_resi_*.root",
        "BcIn": f"{run_tag}_BcIn_resi_*.root",
        "D5": f"{run_tag}_D5_*.root",
    }
    cands = sorted(decode_dir.glob(patterns[kind]))
    return cands[0] if cands else None


def move_misc_to_scratch(dry_run: bool = False):
    data_dir = config.DATA_DIR
    scratch = config.SCRATCH_DIR
    moved = 0

    for entry in sorted(data_dir.iterdir()):
        if not entry.is_file() or entry.is_symlink():
            continue
        dest_dir = None
        if IrregularTpcRe.match(entry.name):
            run_tag = f"run{int(IrregularTpcRe.match(entry.name).group(1)):05d}"
            dest_dir = scratch / run_tag
        elif MiscRunRe.match(entry.name) and not RunLinkRe.match(entry.name):
            run_tag = MiscRunRe.match(entry.name).group(1)
            dest_dir = scratch / run_tag
        if dest_dir is None:
            continue

        dest = dest_dir / entry.name
        print(colored(f"[SCRATCH] {entry.name} -> {dest}", "cyan"))
        if dry_run:
            moved += 1
            continue
        dest_dir.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            dest.unlink()
        shutil.move(str(entry), str(dest))
        moved += 1

    print(colored(f"[SCRATCH] moved {moved} file(s)", "green"))
    return moved


def fix_calibration_symlinks(dry_run: bool = False):
    data_dir = config.DATA_DIR
    fixed = 0

    for link in sorted(data_dir.glob("run*_*.root")):
        if not link.is_symlink():
            continue
        m = RunLinkRe.match(link.name)
        if not m:
            continue
        run_tag, kind = m.group(1), m.group(2)
        final = find_calibration_final(run_tag, kind)
        if final is None:
            continue
        if link.resolve() == final.resolve():
            continue
        if "scratch_root" not in str(link.resolve()) and str(config.DECODE_DIR) in str(link.resolve()):
            continue
        print(colored(f"[LINK] {link.name} -> {final}", "cyan"))
        if not dry_run:
            link.unlink()
            os.symlink(final.resolve(), link)
        fixed += 1

    print(colored(f"[LINK] fixed {fixed} symlink(s)", "green"))
    return fixed


def organize_storage(dry_run: bool = False):
    print(colored("=== Organize group storage ===", "cyan"))
    move_misc_to_scratch(dry_run=dry_run)
    fix_calibration_symlinks(dry_run=dry_run)
    migrate_all_tpchit_roots(run_nums=None, dry_run=dry_run)


def tpchit_era_conf_name(run_num: int) -> str:
    if run_num < 2568:
        return "analyzer_e72_example1.conf"
    if run_num < 3642:
        return "analyzer_e72_example2.conf"
    return "analyzer_e72_example3.conf"


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


def ensure_tpchit_conf(run_num: int) -> str:
    conf_dir = config.PARAM_DIR / "conf"
    sub_conf_dir = conf_dir / config.SUB_DIR
    sub_conf_dir.mkdir(parents=True, exist_ok=True)

    target = sub_conf_dir / f"analyzer_run{run_num:05d}_tpchit.conf"
    tpl_path = conf_dir / TPCHIT_CONF_TEMPLATE
    era_path = conf_dir / tpchit_era_conf_name(run_num)
    digit, cmap = read_digit_cmap(era_path)

    out_lines = []
    for line in tpl_path.read_text().splitlines(keepends=True):
        stripped = line.lstrip()
        if stripped.startswith("DIGIT:"):
            prefix = line[: len(line) - len(stripped)]
            newline = "\n" if line.endswith("\n") else ""
            out_lines.append(f"{prefix}DIGIT:\t\t{digit}{newline}")
        elif stripped.startswith("CMAP:"):
            prefix = line[: len(line) - len(stripped)]
            newline = "\n" if line.endswith("\n") else ""
            out_lines.append(f"{prefix}CMAP:\t\t{cmap}{newline}")
        else:
            out_lines.append(line)

    target.write_text("".join(out_lines))
    return f"param/conf/{config.SUB_DIR}/{target.name}"


def write_tpchit_runlist(path: Path, runs_info: list):
    runlist_dir = config.ANALYZER_DIR / "runmanager/runlist"
    tpl_path = runlist_dir / "tpchit_example.yml"
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(tpl_path) as f_tpl:
        tpl_lines = f_tpl.readlines()

    with open(path, "w") as f_out:
        for line in tpl_lines:
            if line.strip() == "RUN:":
                break
            f_out.write(line)
        f_out.write("RUN:\n")
        for run_info in runs_info:
            f_out.write(f"  {run_info['label']}\n")
            f_out.write(f"    bin: {run_info['bin']}\n")
            f_out.write(f"    conf: {run_info['conf']}\n")
            f_out.write(f"    data: {run_info['data']}\n")
            f_out.write(f"    root: {run_info['root']}\n")
            f_out.write(f"    unit: {run_info['unit']}\n")


def tpchit_batch_runlist_name(run_nums: list) -> str:
    run_nums = sorted(run_nums)
    if len(run_nums) == 1:
        return f"tpchit_run{run_nums[0]:05d}.yml"
    first, last = run_nums[0], run_nums[-1]
    return f"tpchit_{first:05d}-{last:05d}_{len(run_nums)}runs.yml"


def setup_tpchit_run(run_num: int):
    run_tag = f"run{run_num:05d}"
    prefix = "TPC"
    mode_label = "tpchit"

    (config.OUTPUT_DIR / "root" / run_tag).mkdir(parents=True, exist_ok=True)
    (config.SCRATCH_DIR / run_tag).mkdir(parents=True, exist_ok=True)
    (config.DECODE_DIR / run_tag).mkdir(parents=True, exist_ok=True)

    root_name = f"{run_tag}_TPC.root"
    actual_root_abs = config.DECODE_DIR / run_tag / root_name

    symlink_path = config.DATA_DIR / f"{run_tag}_TPC.root"
    ensure_data_symlink(symlink_path, actual_root_abs)

    runlist_dir = config.ANALYZER_DIR / "runmanager/runlist"
    runlist_target_file = runlist_dir / config.SUB_DIR / f"{mode_label}_{run_tag}.yml"
    conf_rel = ensure_tpchit_conf(run_num)

    run_info = {
        "label": f"{run_tag}_tpchit:",
        "bin": "./bin/TPCHit",
        "conf": conf_rel,
        "data": f"rawdata/{run_tag}.dat",
        "root": str(actual_root_abs),
        "unit": 500,
    }

    write_tpchit_runlist(runlist_target_file, [run_info])

    return runlist_target_file, symlink_path, actual_root_abs, run_info


def param_abs_path(rel_path: str) -> Path:
    if rel_path.startswith("param/"):
        return config.PARAM_DIR / rel_path[len("param/"):]
    return config.PARAM_DIR / rel_path


def read_conf_param_path(conf_path: Path, p_key: str):
    for line in conf_path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == p_key:
            return parts[1]
    return None


def find_ref_conf(conf_dir: Path, sub_dir: str, param_run: int, mode_label: str, suffix_head: str):
    primary = conf_dir / sub_dir / f"analyzer_run{param_run:05d}_{mode_label}_{suffix_head}.conf"
    if primary.exists():
        return primary, False
    pattern = f"analyzer_run{param_run:05d}_{mode_label}_*.conf"
    cands = sorted((conf_dir / sub_dir).glob(pattern))
    if cands:
        return cands[0], True
    return None, False


def run_calibration_mode(args):
    SUB_DIR = config.SUB_DIR
    param_run = args.ref if args.ref is not None else args.run_num
    use_ref = args.ref is not None
    run_num = args.run_num

    (config.OUTPUT_DIR / "root" / f"run{run_num:05d}").mkdir(parents=True, exist_ok=True)
    (config.SCRATCH_DIR / f"run{run_num:05d}").mkdir(parents=True, exist_ok=True)
    (config.DECODE_DIR / f"run{run_num:05d}").mkdir(parents=True, exist_ok=True)

    for param_type in ["conf", "USER", "HDPRM", "HDPHC", "DCTDC", "DCDRFT", "DCGEO", "TM"]:
        (config.PARAM_DIR / param_type / SUB_DIR).mkdir(parents=True, exist_ok=True)

    prefix = "Hodo"
    mode_label = "hodo"
    if args.bcout:
        prefix = "BcOut"
        mode_label = "bcout"
    elif args.bcin:
        prefix = "BcIn"
        mode_label = "bcin"
    elif args.d5:
        prefix = "D5"
        mode_label = "d5"

    set1 = set(args.suffix)
    set2 = {"0", "Pi_hdprm", "K_hdprm", "Pi_t0", "K_t0", "Pi_hdphc", "K_hdphc"}
    set3 = {"0", "Pi_tdc", "K_tdc", "Pi_drift", "K_drift", "Pi_resi", "K_resi"}
    set_d5 = {"Pi", "K"}
    FINAL_SUFFIXES = ["hdphc", "resi"]

    if mode_label == "hodo":
        diff = set1 - set2
        if diff or not set1:
            print(f"Error: Invalid suffixes {diff}. Use: 0, Pi/K_hdprm, Pi/K_t0, Pi/K_hdphc")
            sys.exit(1)
    elif mode_label == "d5":
        diff = set1 - set_d5
        if diff or not set1:
            print(f"Error: Invalid suffixes {diff}. Use: Pi, K")
            sys.exit(1)
    else:
        diff = set1 - set3
        if diff or not set1:
            print(f"Error: Invalid suffixes {diff}. Use: 0, Pi/K_tdc, Pi/K_drift, Pi/K_resi")
            sys.exit(1)

    PARAM_DEFS = {
        "USER:":   {"dir": "USER",   "prefix": "UserParam_run",        "tpl": "UserParam_e72_20251104"},
        "HDPRM:":  {"dir": "HDPRM",  "prefix": "HodoParam_run",        "tpl": "HodoParam_e72_example"},
        "HDPHC:":  {"dir": "HDPHC",  "prefix": "HodoPHCParam_run",     "tpl": "HodoPHCParam_e72_example"},
        "DCTDC:":  {"dir": "DCTDC",  "prefix": "DCTdcParam_run",       "tpl": "DCTdcParam_e72_example"},
        "DCDRFT:": {"dir": "DCDRFT", "prefix": "DCDriftParam_run",     "tpl": "DCDriftParam_e72_example.root"},
        "DCGEO:":  {"dir": "DCGEO",  "prefix": "DCGeomParam_run",      "tpl": "DCGeomParam_e72_example"},
        "D5MTX:":  {"dir": "TM",     "prefix": "D5TransferMatrix_run", "tpl": "D5TransferMatrix_example.param"},
    }
    NOSUFFIX_KEYS = {"USER:", "D5MTX:"}
    conf_dir = config.PARAM_DIR / "conf"

    def resolve_param_rel_path(p_key, suffix_head):
        p_info = PARAM_DEFS[p_key]
        run_for_file = run_num if p_key == "D5MTX:" else param_run
        filename = f"{p_info['prefix']}{run_for_file:05d}"
        if p_key not in NOSUFFIX_KEYS:
            filename += f"_{suffix_head}"
        if p_key == "DCDRFT:":
            filename += ".root"

        tier1_rel = f"param/{p_info['dir']}/{SUB_DIR}/{filename}"
        tier1_abs = config.PARAM_DIR / p_info['dir'] / SUB_DIR / filename
        if tier1_abs.exists():
            return tier1_rel

        if p_key == "D5MTX:":
            return f"param/{p_info['dir']}/{p_info['tpl']}"

        tier3_rel = f"param/{p_info['dir']}/{p_info['tpl']}"

        if use_ref:
            ref_conf, ref_conf_fallback = find_ref_conf(
                conf_dir, SUB_DIR, param_run, mode_label, suffix_head,
            )
            if ref_conf is not None:
                borrowed = read_conf_param_path(ref_conf, p_key)
                if borrowed is not None:
                    borrowed_abs = param_abs_path(borrowed)
                    if borrowed_abs.exists():
                        src_note = ref_conf.name
                        if ref_conf_fallback:
                            src_note += " (fallback conf)"
                        print(colored(
                            f"Warning: --ref {param_run} {filename} not found; "
                            f"using {borrowed} from {src_note}",
                            "yellow",
                        ))
                        return borrowed

            tier3_abs = param_abs_path(tier3_rel)
            if tier3_abs.exists():
                print(colored(
                    f"Warning: --ref {param_run} {filename} not found; "
                    f"using template {tier3_rel}",
                    "yellow",
                ))
                return tier3_rel
            print(colored(f"Error: param not found: {tier3_abs}", "red"))
            sys.exit(1)

        return tier3_rel

    d5_mtx_dest = None
    if mode_label == "d5":
        d5_mtx_src = config.PARAM_DIR / "TM" / "D5TransferMatrix_example.param"
        d5_mtx_dest = config.PARAM_DIR / "TM" / SUB_DIR / f"D5TransferMatrix_run{run_num:05d}"
        if not d5_mtx_dest.exists():
            if not d5_mtx_src.exists():
                print(colored(f"Error: TM template not found: {d5_mtx_src}", "red"))
                sys.exit(1)
            shutil.copy(d5_mtx_src, d5_mtx_dest)
        if args.mom is not None:
            set_central_momentum(d5_mtx_dest, args.mom)

    all_runs_info = []
    symlink_path = None

    for suffix in args.suffix:
        suffix_head = suffix.split("_")[0]
        suffix_type = suffix.split("_")[1] if "_" in suffix else suffix
        is_final = suffix_type in FINAL_SUFFIXES or mode_label == "d5"

        conf_target_file = conf_dir / SUB_DIR / f"analyzer_run{run_num:0=5}_{mode_label}_{suffix_head}.conf"
        if not conf_target_file.exists():
            if run_num < 2568:
                example_base = "analyzer_e72_example1.conf"
            elif run_num < 3642:
                example_base = "analyzer_e72_example2.conf"
            else:
                example_base = "analyzer_e72_example3.conf"
            shutil.copy(conf_dir / example_base, conf_target_file)

        buf = []
        has_d5mtx = False
        with open(conf_target_file) as f:
            for line in f:
                s_list = line.split()
                if len(s_list) > 1 and s_list[0] in PARAM_DEFS:
                    p_key = s_list[0]
                    if p_key == "D5MTX:":
                        has_d5mtx = True
                    s_list[1] = resolve_param_rel_path(p_key, suffix_head)
                buf.append(s_list)

        if mode_label == "d5" and not has_d5mtx:
            buf.append(["D5MTX:", resolve_param_rel_path("D5MTX:", suffix_head)])

        with open(conf_target_file, mode='w') as f:
            for l in buf:
                f.write(("\t".join(str(item) for item in l) if l else "") + "\n")

        if suffix_head == suffix_type:
            root_name = f"run{run_num:05d}_{prefix}_{suffix_head}.root"
        else:
            root_name = f"run{run_num:05d}_{prefix}_{suffix_type}_{suffix_head}.root"

        out_dir_abs = (config.DECODE_DIR if is_final else config.SCRATCH_DIR) / f"run{run_num:05d}"
        actual_root_abs = out_dir_abs / root_name

        symlink_name = f"run{run_num:05d}_{prefix}.root"
        symlink_path = config.DATA_DIR / symlink_name
        if symlink_path.exists() or symlink_path.is_symlink():
            symlink_path.unlink()
        try:
            os.symlink(actual_root_abs.resolve(), symlink_path)
        except Exception as e:
            print(f"Warning: Failed to create symlink: {e}")

        option = ""
        if mode_label == "hodo":
            binary = "./bin/Hodoscope"
            unit = 100000
        elif mode_label == "d5":
            binary = "./bin/D5Tracking"
            unit = 20000
            option = "-n 2"
        elif args.bcin:
            binary = "./bin/BcInTracking"
            unit = 50000
            option = "-n 2"
        else:
            binary = "./bin/BcOutTracking"
            unit = 50000
            option = "-n 2"
        all_runs_info.append({
            "label": f"run{run_num:05d}_{suffix}_{mode_label}:",
            "bin": binary,
            "conf": f"param/conf/{SUB_DIR}/{conf_target_file.name}",
            "data": f"rawdata/run{run_num:05d}.dat",
            "root": str(actual_root_abs),
            "unit": unit,
            "option": option,
        })

    runlist_dir = config.ANALYZER_DIR / "runmanager/runlist"
    runlist_target_file = runlist_dir / SUB_DIR / f"{mode_label}_run{run_num:05d}.yml"

    with open(runlist_dir / "myexample.yml") as f_tpl:
        tpl_lines = f_tpl.readlines()

    with open(runlist_target_file, "w") as f_out:
        for line in tpl_lines:
            if line.strip() == "RUN:":
                break
            f_out.write(line)
        f_out.write("RUN:\n")
        for r in all_runs_info:
            f_out.write(f"  {r['label']}\n")
            f_out.write(f"    bin: {r['bin']}\n")
            f_out.write(f"    conf: {r['conf']}\n")
            f_out.write(f"    data: {r['data']}\n")
            f_out.write(f"    root: {r['root']}\n")
            f_out.write(f"    unit: {r['unit']}\n")
            f_out.write(f"    option: {r['option']}\n")

    print(colored(f"\n[SUCCESS] Setup for Run {run_num} {args.suffix} ({prefix})", "green"))
    print(f"  - Params from: run{param_run:05d}" + (" (--ref)" if use_ref else ""))
    if d5_mtx_dest is not None:
        print(f"  - D5MTX: {d5_mtx_dest}")
        if args.mom is not None:
            print(f"  - CentralMomentum: {args.mom} GeV/c")
    print(f"  - Stable Link: {symlink_path}")

    print(colored("\n--- Config File Content ---", "cyan"))
    sample_suffix_head = args.suffix[0].split("_")[0]
    conf_file = config.PARAM_DIR / "conf" / SUB_DIR / f"analyzer_run{run_num:0=5}_{mode_label}_{sample_suffix_head}.conf"
    print(colored(f"File: {conf_file.relative_to(config.ANALYZER_DIR)}", "yellow"))
    with open(conf_file) as f:
        print(f.read())

    print(colored("--- Runlist File Content ---", "cyan"))
    print(colored(f"File: {runlist_target_file.relative_to(config.ANALYZER_DIR)}", "yellow"))
    with open(runlist_target_file) as f:
        print(f.read())
    print("-" * 60)


def main():
    parser = argparse.ArgumentParser(
        prog="make_runlist",
        description="Tool to create analyzer .conf and .yml runlist files with smart storage management.",
    )
    parser.add_argument(
        "run_num", type=int, nargs="?",
        help="Run number (optional with --tpchit --migrate only)",
    )
    parser.add_argument(
        "extra", type=str, nargs="*",
        help="Suffixes (calibration) or additional run numbers (--tpchit)",
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--bcout', action="store_true", help='Set mode to BcOut calibration')
    group.add_argument('--bcin', action="store_true", help='Set mode to BcIn calibration')
    group.add_argument('--d5', action="store_true", help='Set mode to D5Tracking')
    group.add_argument('--tpchit', action="store_true", help='Set mode to TPCHit decode')
    parser.add_argument("--ref", type=int, default=None, help="Use calibrated params from this run")
    parser.add_argument("--mom", type=float, default=None, help="D5 central momentum in GeV/c")
    parser.add_argument(
        "--migrate", action="store_true",
        help="Migrate existing tpc_run*.root files to root/runNNNNN/runNNNNN_TPC.root",
    )
    parser.add_argument(
        "--organize", action="store_true",
        help="Move irregular files to scratch and fix calibration symlinks",
    )
    parser.add_argument("--dry-run", action="store_true", help="Show actions without moving files")

    args = parser.parse_args()

    if args.organize:
        organize_storage(dry_run=args.dry_run)
        return

    if args.tpchit:
        run_nums = []
        if args.run_num is not None:
            run_nums.append(args.run_num)
        for x in args.extra:
            if x.isdigit():
                run_nums.append(int(x))
            else:
                print(colored(f"Error: unknown argument '{x}' for --tpchit", "red"))
                sys.exit(1)
        suffix = []
        if not run_nums and not args.migrate:
            print(colored("Error: --tpchit requires at least one run number", "red"))
            sys.exit(1)
    else:
        if args.run_num is None:
            print(colored("Error: run number is required", "red"))
            sys.exit(1)
        run_nums = [args.run_num]
        suffix = args.extra

    if args.mom is not None and not args.d5:
        print(colored("Error: --mom is only valid with --d5", "red"))
        sys.exit(1)

    if args.tpchit:
        check_cobo_consistency()
        if args.migrate:
            migrate_all_tpchit_roots(run_nums or None, dry_run=args.dry_run)
            if args.dry_run:
                return
        all_run_infos = []
        for run_num in sorted(set(run_nums)):
            yml, link, root, run_info = setup_tpchit_run(run_num)
            all_run_infos.append(run_info)
            era = tpchit_era_conf_name(run_num)
            print(colored(f"\n[SUCCESS] TPCHit setup for run {run_num:05d}", "green"))
            print(f"  - Conf: {run_info['conf']} (DIGIT/CMAP from {era})")
            print(f"  - Root: {root}")
            print(f"  - Link: {link}")
            print(colored(f"  - Runlist: {yml.relative_to(config.ANALYZER_DIR)}", "cyan"))
            print(colored("\n--- Runlist File Content ---", "cyan"))
            with open(yml) as f:
                print(f.read())

        if len(all_run_infos) > 1:
            runlist_dir = config.ANALYZER_DIR / "runmanager/runlist"
            batch_name = tpchit_batch_runlist_name(sorted(set(run_nums)))
            batch_yml = runlist_dir / TPCHIT_SCRATCH_RUNLIST_DIR / batch_name
            write_tpchit_runlist(batch_yml, all_run_infos)
            print(colored(
                f"\n[BATCH] Combined runlist for parallel decode ({len(all_run_infos)} runs)",
                "green",
            ))
            print(colored(
                f"  - {batch_yml.relative_to(config.ANALYZER_DIR)}",
                "cyan",
            ))
            print(colored("\n--- Batch Runlist File Content ---", "cyan"))
            with open(batch_yml) as f:
                print(f.read())
        return

    if args.migrate:
        print(colored("Error: --migrate requires --tpchit", "red"))
        sys.exit(1)

    class CalibArgs:
        pass

    calib = CalibArgs()
    calib.run_num = run_nums[0]
    calib.suffix = suffix
    calib.bcout = args.bcout
    calib.bcin = args.bcin
    calib.d5 = args.d5
    calib.ref = args.ref
    calib.mom = args.mom
    run_calibration_mode(calib)


if __name__ == "__main__":
    main()
