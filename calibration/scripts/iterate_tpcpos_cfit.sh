#!/usr/bin/env bash
# TPCPOS iteration from calibration data only (beam-through + multi-pion):
#   BT17 Helix(map_k) (full via myrun, or 1/3 skim via btskim/bt_chain.sh with --bt-skim)
#   + HTOF6 skim Helix(map_k)
#   -> tpcpos_constrained_fit (BT radius-constrained or, with free, inclusive residual; HTOF6 inclusive
#      residual fills pads without a BT update) -> map_{k+1} = map_k + update -> ...
# K-p / Lambda checks (lamk_chain.sh / lam_chain.sh) only benchmark each map; they do not
# feed back into the fit.
# Spec: docs/tpc/TPCPOS_ITER_PREP.ja.md
#
# Usage (blocks; run under nohup/tmux):
#   iterate_tpcpos_cfit.sh --start LABEL --n N [--start-bt-hadd FILE] [--alpha A]
#                          [--ref-mode self|pair|soft|free] [--sigma-p S] [--ref-weight W] [--prefix s] [--first K]
#                          [--no-bench] [--bt-skim] [--zero-layers L0-L1]
#     LABEL : map label; map at bt_recon/maps/LABEL.map, BT17 Helix in bt_recon/LABEL
#             (or --start-bt-hadd), HTOF6 Helix by htof/htof_chain.sh LABEL.
#   New maps are <prefix>K, <prefix>K+1, ... (K = --first, default 1).
#   Python runs from the real conda env (~/.conda/envs/work), not the node-local envfs copy.
#
# Shared e72 map / standard ROOT files are not touched.

set -euo pipefail

START=""
START_BT_HADD=""
N_ITER=2
ALPHA=0.7
REF_MODE=self
REF_W=1.0
SIGMA_P=0.01
PREFIX="s"
FIRST=1
BT_SKIM=0
ZERO_LAYERS=""
BENCH=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --start) START="$2"; shift 2 ;;
    --start-bt-hadd) START_BT_HADD="$2"; shift 2 ;;
    --n) N_ITER="$2"; shift 2 ;;
    --alpha) ALPHA="$2"; shift 2 ;;
    --ref-mode) REF_MODE="$2"; shift 2 ;;
    --ref-weight) REF_W="$2"; shift 2 ;;
    --sigma-p) SIGMA_P="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --first) FIRST="$2"; shift 2 ;;
    --no-bench) BENCH=0; shift ;;
    --bt-skim) BT_SKIM=1; shift ;;
    --zero-layers) ZERO_LAYERS="$2"; shift 2 ;;
    *) echo "[Error] unknown arg: $1" >&2; exit 1 ;;
  esac
done
[[ -n "${START}" ]] || { echo "[Error] --start LABEL required" >&2; exit 1; }

A=/gpfs/home/had/sryuta/analyzer/JPARC2025E72
SCR=${A}/myanalysis/analysis/scripts
TP=/group/had/sks/Users/sryuta/tmp/debug/tpcpos
T=${TP}/bt_recon
H=${TP}/hadd
L=${TP}/lamk
CF=${TP}/cfit
HT=${TP}/htof
BS=${TP}/btskim

export PATH="${HOME}/.conda/envs/work/bin:${PATH}"

log() { echo "[$(date '+%m-%d %T')] $*"; }

wait_bt() {
  local label="$1"
  if [[ ${BT_SKIM} -eq 1 ]]; then
    until grep -q "bt done\|failures" ${BS}/chain_${label}.log 2>/dev/null; do sleep 120; done
    grep -q "bt done" ${BS}/chain_${label}.log || { log "[Error] BT17 skim helix failed (${label})"; exit 1; }
    ln -sfn ${BS}/bt17s_${label}.root ${H}/bt17_${label}.root
    return 0
  fi
  while pgrep -u "$(id -u)" -f "myrun.py ${T}/dst_tpchelix_bt17_${label}.yml" > /dev/null; do sleep 120; done
  local n
  n=$(ls ${T}/${label}/run*_TPCHelix.root 2>/dev/null | wc -l)
  [[ ${n} -ge 17 ]] || { log "[Error] ${label}: only ${n}/17 BT Helix files"; exit 1; }
  hadd -f ${H}/bt17_${label}.root ${T}/${label}/run*_TPCHelix.root > ${H}/bt17_${label}_hadd.log 2>&1
}

wait_htof() {
  local label="$1"
  until grep -q "htof done\|failures" ${HT}/chain_${label}.log 2>/dev/null; do sleep 120; done
  grep -q "htof done" ${HT}/chain_${label}.log || { log "[Error] HTOF6 helix failed (${label})"; exit 1; }
}

start_next() {
  local label="$1"
  if [[ ${BT_SKIM} -eq 1 ]]; then
    (cd ${BS} && nohup ./bt_chain.sh ${label} ${CF}/${label}.map > ${BS}/chain_${label}.log 2>&1 &)
  else
    ${T}/make_cand.sh "${label}" ${CF}/${label}.map > /dev/null
    (cd ${A} && nohup ./runmanager/myrun.py ${T}/dst_tpchelix_bt17_${label}.yml \
       > ${T}/myrun_${label}.log 2>&1 < /dev/null &)
  fi
  (cd ${HT} && nohup ./htof_chain.sh ${label} ${CF}/${label}.map > ${HT}/chain_${label}.log 2>&1 &)
  if [[ ${BENCH} -eq 1 ]]; then
    (cd ${L} && nohup ./lamk_chain.sh ${label} ${CF}/${label}.map > ${L}/chain_${label}.log 2>&1 &)
    (cd ${L} && RUN=03778 nohup ./lamk_chain.sh ${label} ${CF}/${label}.map > ${L}/chain3778_${label}.log 2>&1 &)
    (cd ${L} && nohup ./lam_chain.sh ${label} ${CF}/${label}.map > ${L}/chainlam_${label}.log 2>&1 &)
  fi
}

qa() {
  local tag="$1"; shift
  local bt=() ht=()
  for lb in "$@"; do
    bt+=("${lb}=${H}/bt17_${lb}.root")
    ht+=(--htof "${lb}=${HT}/htof6_${lb}.root")
  done
  python3 ${SCR}/tpcpos_iter_eval.py "${bt[@]}" "${ht[@]}" --tag ${tag} \
    > ${TP}/cfit_iter_eval_${tag}.log 2>&1 || log "[WARN] iter_eval failed"
  sed -n '/^fit/,/^charge-pair/p' ${TP}/cfit_iter_eval_${tag}.log
}

cur="${START}"
labels=("${START}")
log "start=${START} n=${N_ITER} alpha=${ALPHA} ref=${REF_MODE} w=${REF_W} sigma_p=${SIGMA_P}"
for ((k = FIRST; k < FIRST + N_ITER; k++)); do
  next="${PREFIX}${k}"
  if [[ ${k} -eq ${FIRST} && -n "${START_BT_HADD}" ]]; then
    ln -sfn "${START_BT_HADD}" ${H}/bt17_${cur}.root
  else
    log "wait BT17 Helix (${cur})"
    wait_bt "${cur}"
  fi
  log "wait HTOF6 Helix (${cur})"
  wait_htof "${cur}"
  parent=${CF}/${cur}.map
  [[ -f ${T}/maps/${cur}.map ]] && parent=${T}/maps/${cur}.map
  python3 ${SCR}/tpcpos_constrained_fit.py ${H}/bt17_${cur}.root \
    --parent ${parent} --out ${CF}/${next}.map --alpha ${ALPHA} \
    --ref-mode ${REF_MODE} --ref-weight ${REF_W} --sigma-p ${SIGMA_P} \
    --htof-inclusive ${HT}/htof6_${cur}.root ${ZERO_LAYERS:+--zero-layers ${ZERO_LAYERS}} \
    > ${CF}/${next}.log 2>&1
  log "fit ${cur} -> ${next}: $(grep -E 'inclusive-filled|tracks=' ${CF}/${next}.log | tr '\n' ' ')"
  start_next "${next}"
  log "submitted BT17 + HTOF6 Helix (and benchmarks) for ${next}"
  log "QA up to ${cur}:"
  qa "${PREFIX}_upto_${cur}" "${labels[@]}"
  labels+=("${next}")
  cur="${next}"
done
log "last map ${cur}: wait Helix for final QA"
wait_bt "${cur}"
wait_htof "${cur}"
qa "${PREFIX}_upto_${cur}" "${labels[@]}"
log "done"
