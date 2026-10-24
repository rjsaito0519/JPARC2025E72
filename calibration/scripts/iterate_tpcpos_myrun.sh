#!/usr/bin/env bash
# TPCPOS complementary iter: Helix → pos per sample (BT then HTOF), not back-to-back pos.
# Spec: docs/tpc/TPCPOS_ITER_PREP.ja.md
#
# Location: myanalysis/calibration/scripts/iterate_tpcpos_myrun.sh
# Usage (inside tmux):
#   cd ~/JPARC2025E72/calibration/scripts
#   ./iterate_tpcpos_myrun.sh --n 2
#   ./iterate_tpcpos_myrun.sh --n 1 --skip-pos     # Helix only (map already OK)
#   ./iterate_tpcpos_myrun.sh --n 1 --bt-only
#   ./iterate_tpcpos_myrun.sh --dry-run
#   ./iterate_tpcpos_myrun.sh --n 1 --skip-helix --parent-bt MAP --parent-htof MAP
#
# pos adds its update onto the map that reconstructed the Helix (--parent).
#
# Per iter (default, complementary):
#   map backup
#   BT:   myrun → hadd → (QA pos --debug) → pos write (no --ave)
#   HTOF: myrun → hadd → (QA pos --debug) → pos write (--ave)
#   → optional Enter gate
#
# Fixed runlists (dst_create_runlist --tpchelix; conf = *_tpcpos.conf → e72 map):
#   runmanager/runlist/scratch/dst_tpchelix_bt17_must.yml
#   runmanager/runlist/scratch/dst_tpchelix_htof2601_6.yml

set -euo pipefail

N_ITER=2
DRY_RUN=0
SKIP_POS=0
SKIP_HELIX=0
BT_ONLY=0
HTOF_ONLY=0
DEBUG_QA=1
NO_GATE=0
THRESHOLD=""
PARENT_BT=""
PARENT_HTOF=""
ANALYZER_DIR="/home/had/sryuta/analyzer/JPARC2025E72"
MYANALYSIS_DIR="${ANALYZER_DIR}/myanalysis"
GROUP_ROOT="${ANALYZER_DIR}/group_root"
TMP_DIR="/group/had/sks/Users/sryuta/tmp/debug/tpcpos"
HADD_DIR="${TMP_DIR}/hadd"
MAP_DIR="${TMP_DIR}/maps"
MAP_REL="param/TPCPOS/e72/TPCPositionCorrectionMap"
YML_BT="runmanager/runlist/scratch/dst_tpchelix_bt17_must.yml"
YML_HTOF="runmanager/runlist/scratch/dst_tpchelix_htof2601_6.yml"
ADD_BT="${GROUP_ROOT}/add_helix_bt17.sh"
ADD_HTOF="${GROUP_ROOT}/add_helix_htof2601.sh"

usage() {
  sed -n '2,24p' "$0" | sed 's/^# \?//'
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --n) N_ITER="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --skip-pos) SKIP_POS=1; shift ;;
    --skip-helix) SKIP_HELIX=1; shift ;;
    --bt-only) BT_ONLY=1; shift ;;
    --htof-only) HTOF_ONLY=1; shift ;;
    --no-qa) DEBUG_QA=0; shift ;;
    --yes|--no-gate) NO_GATE=1; shift ;;
    --threshold|--min-entries) THRESHOLD="$2"; shift 2 ;;
    --parent-bt) PARENT_BT="$2"; shift 2 ;;
    --parent-htof) PARENT_HTOF="$2"; shift 2 ;;
    --analyzer-dir) ANALYZER_DIR="$2"; MYANALYSIS_DIR="${ANALYZER_DIR}/myanalysis"; GROUP_ROOT="${ANALYZER_DIR}/group_root"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "[Error] unknown arg: $1" >&2; usage ;;
  esac
done

if [[ "${BT_ONLY}" -eq 1 && "${HTOF_ONLY}" -eq 1 ]]; then
  echo "[Error] use only one of --bt-only / --htof-only" >&2
  exit 1
fi

RUN_TPC="${MYANALYSIS_DIR}/run_tpc.py"
MYRUN="${ANALYZER_DIR}/runmanager/myrun.py"
MAP_ABS="${ANALYZER_DIR}/${MAP_REL}"
YML_BT_ABS="${ANALYZER_DIR}/${YML_BT}"
YML_HTOF_ABS="${ANALYZER_DIR}/${YML_HTOF}"

mkdir -p "${HADD_DIR}" "${MAP_DIR}"

require_file() {
  local p="$1"
  if [[ ! -f "${p}" ]]; then
    echo "[Error] missing: ${p}" >&2
    exit 1
  fi
}

require_file "${RUN_TPC}"
require_file "${MYRUN}"
require_file "${MAP_ABS}"
require_file "${YML_BT_ABS}"
require_file "${YML_HTOF_ABS}"
require_file "${ADD_BT}"
require_file "${ADD_HTOF}"

check_yml() {
  local yml="$1"
  local tag="$2"
  if ! grep -q 'unit: 2000' "${yml}"; then
    echo "[Error] ${tag} yml missing unit: 2000: ${yml}" >&2
    exit 1
  fi
  if ! grep -q '_tpcpos.conf' "${yml}"; then
    echo "[Error] ${tag} yml should use *_tpcpos.conf (e72 map): ${yml}" >&2
    exit 1
  fi
}
check_yml "${YML_BT_ABS}" "BT"
check_yml "${YML_HTOF_ABS}" "HTOF"

run_cmd() {
  echo "[cmd] $*"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    return 0
  fi
  "$@"
}

LAST_BACKUP=""
backup_map() {
  local tag="$1"
  local stamp
  stamp="$(date +%Y%m%d_%H%M%S)"
  local dest="${MAP_DIR}/TPCPositionCorrectionMap_e72_${tag}_${stamp}"
  echo "[map] backup → ${dest}"
  LAST_BACKUP="${dest}"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    return 0
  fi
  cp -a "${MAP_ABS}" "${dest}"
}

# pos adds onto --parent = the map snapshot taken right before that Helix.
run_pos() {
  local root="$1"
  local ave_flag="${2:-}"
  local debug_flag="${3:-}"
  local parent="${4:-}"
  local cmd=("${RUN_TPC}" "${root}" pos)
  if [[ -n "${parent}" ]]; then
    cmd+=(--parent "${parent}")
  fi
  if [[ -n "${ave_flag}" ]]; then
    cmd+=(--ave)
  fi
  if [[ -n "${debug_flag}" ]]; then
    cmd+=(--debug)
  fi
  if [[ -n "${THRESHOLD}" ]]; then
    cmd+=(--threshold "${THRESHOLD}")
  fi
  run_cmd "${cmd[@]}"
}

run_myrun() {
  local yml="$1"
  local tag="$2"
  echo "[myrun] ${tag}: ${yml}"
  echo "  (blocking — wait until jobs finish; optional: monitor.py in another pane)"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[dry-run] would: cd ${ANALYZER_DIR} && ./runmanager/myrun.py ${yml#"${ANALYZER_DIR}/"}"
    return 0
  fi
  (
    cd "${ANALYZER_DIR}"
    ./runmanager/myrun.py "${yml#"${ANALYZER_DIR}/"}"
  )
}

do_hadd_bt() {
  local out="$1"
  echo "[hadd] BT17 → ${out}"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[dry-run] OUT=${out} ${ADD_BT}"
    return 0
  fi
  OUT="${out}" bash "${ADD_BT}"
}

do_hadd_htof() {
  local out="$1"
  echo "[hadd] HTOF6 → ${out}"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[dry-run] OUT=${out} ${ADD_HTOF}"
    return 0
  fi
  OUT="${out}" bash "${ADD_HTOF}"
}

resolve_prev_hadd() {
  local kind="$1"
  local prev="$2"
  local cand
  for cand in \
    "${HADD_DIR}/${kind}_iter${prev}.root" \
    "${HADD_DIR}/${kind}_iter0.root" \
    "${GROUP_ROOT}/${kind}_tpcpos_helix.root"
  do
    if [[ -f "${cand}" ]]; then
      echo "${cand}"
      return 0
    fi
  done
  return 1
}

# Helix (current map) → hadd → optional QA → optional map-updating pos.
leg_bt() {
  local i="$1"
  local prev=$((i - 1))
  local bt_out="${HADD_DIR}/bt17_iter${i}.root"
  local bt_for_pos=""

  local bt_parent="${PARENT_BT}"

  echo "---- BT leg (iter ${i}) ----"
  if [[ "${SKIP_HELIX}" -eq 0 ]]; then
    backup_map "iter${i}_bt_helix_parent"
    bt_parent="${LAST_BACKUP}"
    run_myrun "${YML_BT_ABS}" "BT17"
    do_hadd_bt "${bt_out}"
    bt_for_pos="${bt_out}"
    if [[ "${DEBUG_QA}" -eq 1 ]]; then
      run_pos "${bt_out}" "" --debug "${bt_parent}"
    fi
  else
    echo "[helix] BT skipped (--skip-helix)"
    bt_for_pos="$(resolve_prev_hadd bt17 "${prev}" || true)"
  fi

  if [[ "${SKIP_POS}" -eq 0 ]]; then
    if [[ -z "${bt_for_pos}" ]]; then
      echo "[Error] no BT hadd for pos (Helix this iter or ${HADD_DIR}/bt17_iter${prev}.root)." >&2
      exit 1
    fi
    echo "[pos] BT write input=${bt_for_pos}"
    if [[ -z "${bt_parent}" ]]; then
      echo "[Error] --skip-helix needs --parent-bt <map used for that Helix>." >&2
      exit 1
    fi
    run_pos "${bt_for_pos}" "" "" "${bt_parent}"
    backup_map "iter${i}_after_bt_pos"
  else
    echo "[pos] BT skipped (--skip-pos)"
  fi
}

leg_htof() {
  local i="$1"
  local prev=$((i - 1))
  local htof_out="${HADD_DIR}/htof2601_iter${i}.root"
  local htof_for_pos=""

  local htof_parent="${PARENT_HTOF}"

  echo "---- HTOF leg (iter ${i}; map should include BT pos if not --skip-pos) ----"
  if [[ "${SKIP_HELIX}" -eq 0 ]]; then
    backup_map "iter${i}_htof_helix_parent"
    htof_parent="${LAST_BACKUP}"
    run_myrun "${YML_HTOF_ABS}" "HTOF6"
    do_hadd_htof "${htof_out}"
    htof_for_pos="${htof_out}"
    if [[ "${DEBUG_QA}" -eq 1 ]]; then
      run_pos "${htof_out}" "" --debug "${htof_parent}"
    fi
  else
    echo "[helix] HTOF skipped (--skip-helix)"
    htof_for_pos="$(resolve_prev_hadd htof2601 "${prev}" || true)"
  fi

  if [[ "${SKIP_POS}" -eq 0 ]]; then
    if [[ -z "${htof_for_pos}" ]]; then
      echo "[Error] no HTOF hadd for pos." >&2
      exit 1
    fi
    echo "[pos] HTOF write --ave input=${htof_for_pos}"
    if [[ -z "${htof_parent}" ]]; then
      echo "[Error] --skip-helix needs --parent-htof <map used for that Helix>." >&2
      exit 1
    fi
    run_pos "${htof_for_pos}" --ave "" "${htof_parent}"
    backup_map "iter${i}_after_htof_pos"
  else
    echo "[pos] HTOF skipped (--skip-pos)"
  fi
}

echo "=== TPCPOS iterate_tpcpos_myrun n=${N_ITER} complementary Helix→pos dry=${DRY_RUN} ==="
echo "  map     : ${MAP_ABS}"
echo "  yml BT  : ${YML_BT}"
echo "  yml HTOF: ${YML_HTOF}"
echo "  hadd dir: ${HADD_DIR}"
echo "  order   : BT(Helix→pos) then HTOF(Helix→pos)  [not back-to-back pos]"
if [[ -z "${TMUX:-}" ]]; then
  echo "[WARN] TMUX not set — prefer: tmux new -s tpcpos_iter"
fi

for ((i = 1; i <= N_ITER; i++)); do
  echo ""
  echo "---------- iter ${i}/${N_ITER} ----------"
  backup_map "iter${i}_pre"

  if [[ "${HTOF_ONLY}" -eq 0 ]]; then
    leg_bt "${i}"
  fi
  if [[ "${BT_ONLY}" -eq 0 ]]; then
    leg_htof "${i}"
  fi

  if [[ "${i}" -lt "${N_ITER}" && "${NO_GATE}" -eq 0 && "${DRY_RUN}" -eq 0 ]]; then
    echo ""
    echo "[gate] iter ${i} done. Review QA / map, then Enter for next (or Ctrl-C to stop)."
    read -r _
  fi
done

echo ""
echo "=== done ${N_ITER} iter(s). map=${MAP_ABS} ==="
echo "  Tip: map already written, Helix only → --skip-pos"
echo "  Monitor (optional): cd ${ANALYZER_DIR} && ./runmanager/monitor.py runmanager/stat/<yml>.json"
