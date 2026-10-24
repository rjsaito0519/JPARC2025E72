#!/usr/bin/env bash
# TPCPOS manual iter: BT must17 hadd → pos (no --ave) → HTOF 2601系 hadd → pos (--ave)
# Shared map: param/TPCPOS/e72/TPCPositionCorrectionMap
# Spec: docs/tpc/TPCPOS_ITER_PREP.ja.md
#
# Location: myanalysis/calibration/scripts/iterate_tpcpos_bt_htof.sh
# (hadd+pos only; prefer iterate_tpcpos_myrun.sh for Helix+myrun)
#
# Usage:
#   iterate_tpcpos_bt_htof.sh --n N [--do-recon] [--debug-pos] [--dry-run]
#                            [--threshold N] [--root-dir DIR] [--analyzer-dir DIR]
#
# Defaults: --skip-recon (do not run helix recon).
# --dry-run: print steps only (no hadd / no pos).

set -euo pipefail

N_ITER=1
DO_RECON=0
DEBUG_POS=0
DRY_RUN=0
THRESHOLD=""
ROOT_DIR="/group/had/sks/Users/sryuta/JPARC2025E72/root"
ANALYZER_DIR="/home/had/sryuta/analyzer/JPARC2025E72"
TMP_DIR="/group/had/sks/Users/sryuta/tmp/debug/tpcpos"
HADD_DIR="${TMP_DIR}/hadd"
MAP_REL="param/TPCPOS/e72/TPCPositionCorrectionMap"

# Sufficiency 2026-10-17: must 17 only (optional BT omitted; 2531 dropped)
BT_RUNS=(2489 2491 2494 2500 2502 2504 2509 2512 2514 2516 2518 2520 2522 2524 2527 2535 2537)
HTOF_RUNS=(2601 2602 2603 2604 2606 2607)

usage() {
  sed -n '2,12p' "$0" | sed 's/^# \?//'
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --n) N_ITER="$2"; shift 2 ;;
    --do-recon) DO_RECON=1; shift ;;
    --skip-recon) DO_RECON=0; shift ;;
    --debug-pos) DEBUG_POS=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    --threshold|--min-entries) THRESHOLD="$2"; shift 2 ;;
    --root-dir) ROOT_DIR="$2"; shift 2 ;;
    --analyzer-dir) ANALYZER_DIR="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "[Error] unknown arg: $1" >&2; usage ;;
  esac
done

RUN_TPC="${ANALYZER_DIR}/myanalysis/run_tpc.py"
MAP_ABS="${ANALYZER_DIR}/${MAP_REL}"
mkdir -p "${HADD_DIR}"

if [[ ! -x "${RUN_TPC}" && ! -f "${RUN_TPC}" ]]; then
  echo "[Error] run_tpc not found: ${RUN_TPC}" >&2
  exit 1
fi
if [[ ! -f "${MAP_ABS}" ]]; then
  echo "[Error] shared map missing: ${MAP_ABS}" >&2
  echo "  Init: mkdir -p param/TPCPOS/e72 && cp param/TPCPOS/TPCPositionCorrectionMap_0 param/TPCPOS/e72/TPCPositionCorrectionMap" >&2
  exit 1
fi

helix_path() {
  local run="$1"
  printf "%s/run%05d/run%05d_TPCHelix.root" "${ROOT_DIR}" "${run}" "${run}"
}

collect_existing() {
  local -n _runs=$1
  local -n _out=$2
  _out=()
  local r p
  for r in "${_runs[@]}"; do
    p="$(helix_path "${r}")"
    if [[ -f "${p}" ]]; then
      _out+=("${p}")
    else
      echo "[WARN] missing Helix, skip run ${r}: ${p}" >&2
    fi
  done
}

do_hadd() {
  local out="$1"
  shift
  local inputs=("$@")
  if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "[Error] no inputs for hadd → ${out}" >&2
    return 1
  fi
  echo "[hadd] ${#inputs[@]} files → ${out}"
  # shellcheck disable=SC2086
  hadd -f "${out}" "${inputs[@]}"
}

run_pos() {
  local root="$1"
  local ave_flag="${2:-}"
  local cmd=("${RUN_TPC}" "${root}" pos)
  if [[ -n "${ave_flag}" ]]; then
    cmd+=(--ave)
  fi
  if [[ "${DEBUG_POS}" -eq 1 ]]; then
    cmd+=(--debug)
  fi
  if [[ -n "${THRESHOLD}" ]]; then
    cmd+=(--threshold "${THRESHOLD}")
  fi
  echo "[pos] ${cmd[*]}"
  "${cmd[@]}"
}

echo "=== TPCPOS iter BT→HTOF n=${N_ITER} do_recon=${DO_RECON} debug_pos=${DEBUG_POS} dry_run=${DRY_RUN} ==="
echo "  shared map: ${MAP_ABS}"
echo "  root dir  : ${ROOT_DIR}"
echo "  BT runs   : ${#BT_RUNS[@]} (must17)"

for ((i = 1; i <= N_ITER; i++)); do
  echo ""
  echo "---------- iter ${i}/${N_ITER} ----------"

  bt_files=()
  collect_existing BT_RUNS bt_files
  bt_hadd="${HADD_DIR}/bt17_iter${i}.root"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[dry-run] would hadd ${#bt_files[@]} BT files → ${bt_hadd}"
    echo "[dry-run] would pos (no --ave) → ${MAP_REL}"
  else
    do_hadd "${bt_hadd}" "${bt_files[@]}"
    run_pos "${bt_hadd}" ""
  fi

  htof_files=()
  collect_existing HTOF_RUNS htof_files
  htof_hadd="${HADD_DIR}/htof2601_iter${i}.root"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[dry-run] would hadd ${#htof_files[@]} HTOF files → ${htof_hadd}"
    echo "[dry-run] would pos --ave → ${MAP_REL}"
  else
    do_hadd "${htof_hadd}" "${htof_files[@]}"
    run_pos "${htof_hadd}" --ave
  fi

  if [[ "${DO_RECON}" -eq 1 ]]; then
    echo "[recon] not wired yet (job-gated). Would smoke run 2489 / 2601 with ${MAP_REL}"
  else
    echo "[recon] skipped (--skip-recon default; pass --do-recon when queue is free)"
  fi
done

echo ""
echo "=== done ${N_ITER} iter(s). map target: ${MAP_ABS} ==="
echo "  Note: use --debug-pos for analysis/QA without map write."
