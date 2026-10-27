#!/usr/bin/env bash
# Iterate the TPC hit-position correction map (TPCPOS) from calibration data only
# (beam-through runs + multi-pion runs):
#   beam-through Helix(map_k) (17 runs; full via myrun, or a 1/3 sample via btskim/bt_chain.sh with --bt-skim)
#   + multi-pion Helix(map_k) (HTOF trigger runs 2601-2607, skim)
#   -> tpc_position_update.py (beam-through residual, radius constrained per --momentum-constraint;
#      the multi-pion inclusive residual fills pads without a beam-through update)
#   -> map_{k+1} = map_k + update -> ...
# K-p / Lambda checks (lamk_chain.sh / lam_chain.sh) only check each map; they do not
# feed back into the fit.
#
# Usage (blocks; run under nohup/tmux):
#   iterate_tpc_position.sh --from LABEL --iterations N --name NAME --momentum-constraint MODE
#                           [--momentum-width S] [--pair-weight W] [--update-fraction F]
#                           [--first-iteration K] [--from-bt-hadd FILE]
#                           [--bt-skim] [--zero-layers L0-L1] [--no-kp-check]
#     LABEL : label of the starting map; map at bt_recon/maps/LABEL.map or cfit/LABEL.map,
#             beam-through Helix in bt_recon/LABEL (or --from-bt-hadd), multi-pion Helix by
#             htof/htof_chain.sh LABEL.
#     MODE  : none | soft | fixed | charge-pair (see tpc_position_update.py).
#   New maps are NAME_iterK, NAME_iterK+1, ... (K = --first-iteration, default 1).
#   Python runs from the real conda env (~/.conda/envs/work), not the node-local envfs copy.
#
# Shared e72 map / standard ROOT files are not touched.

set -euo pipefail

START=""
START_BT_HADD=""
N_ITER=""
NAME=""
MOM_CONSTRAINT=""
MOM_WIDTH=0.01
PAIR_W=0.5
FRACTION=0.7
FIRST=1
BT_SKIM=0
ZERO_LAYERS=""
BENCH=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) START="$2"; shift 2 ;;
    --from-bt-hadd) START_BT_HADD="$2"; shift 2 ;;
    --iterations) N_ITER="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --momentum-constraint) MOM_CONSTRAINT="$2"; shift 2 ;;
    --momentum-width) MOM_WIDTH="$2"; shift 2 ;;
    --pair-weight) PAIR_W="$2"; shift 2 ;;
    --update-fraction) FRACTION="$2"; shift 2 ;;
    --first-iteration) FIRST="$2"; shift 2 ;;
    --no-kp-check) BENCH=0; shift ;;
    --bt-skim) BT_SKIM=1; shift ;;
    --zero-layers) ZERO_LAYERS="$2"; shift 2 ;;
    *) echo "[Error] unknown arg: $1" >&2; exit 1 ;;
  esac
done
for req in START:--from N_ITER:--iterations NAME:--name MOM_CONSTRAINT:--momentum-constraint; do
  var=${req%%:*}
  [[ -n "${!var}" ]] || { echo "[Error] ${req#*:} is required" >&2; exit 1; }
done
case "${MOM_CONSTRAINT}" in
  none|soft|fixed|charge-pair) ;;
  *) echo "[Error] --momentum-constraint must be none, soft, fixed or charge-pair" >&2; exit 1 ;;
esac

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
    grep -q "bt done" ${BS}/chain_${label}.log || { log "[Error] beam-through skim Helix failed (${label})"; exit 1; }
    ln -sfn ${BS}/bt17s_${label}.root ${H}/bt17_${label}.root
    return 0
  fi
  while pgrep -u "$(id -u)" -f "myrun.py ${T}/dst_tpchelix_bt17_${label}.yml" > /dev/null; do sleep 120; done
  local n
  n=$(ls ${T}/${label}/run*_TPCHelix.root 2>/dev/null | wc -l)
  [[ ${n} -ge 17 ]] || { log "[Error] ${label}: only ${n}/17 beam-through Helix files"; exit 1; }
  hadd -f ${H}/bt17_${label}.root ${T}/${label}/run*_TPCHelix.root > ${H}/bt17_${label}_hadd.log 2>&1
}

wait_htof() {
  local label="$1"
  until grep -q "htof done\|failures" ${HT}/chain_${label}.log 2>/dev/null; do sleep 120; done
  grep -q "htof done" ${HT}/chain_${label}.log || { log "[Error] multi-pion Helix failed (${label})"; exit 1; }
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
log "from=${START} iterations=${N_ITER} name=${NAME} momentum-constraint=${MOM_CONSTRAINT} width=${MOM_WIDTH} pair-weight=${PAIR_W} update-fraction=${FRACTION}"
for ((k = FIRST; k < FIRST + N_ITER; k++)); do
  next="${NAME}_iter${k}"
  if [[ ${k} -eq ${FIRST} && -n "${START_BT_HADD}" ]]; then
    ln -sfn "${START_BT_HADD}" ${H}/bt17_${cur}.root
  else
    log "wait beam-through Helix (${cur})"
    wait_bt "${cur}"
  fi
  log "wait multi-pion Helix (${cur})"
  wait_htof "${cur}"
  parent=${CF}/${cur}.map
  [[ -f ${T}/maps/${cur}.map ]] && parent=${T}/maps/${cur}.map
  python3 ${SCR}/tpc_position_update.py ${H}/bt17_${cur}.root \
    --parent ${parent} --out ${CF}/${next}.map --update-fraction ${FRACTION} \
    --momentum-constraint ${MOM_CONSTRAINT} --momentum-width ${MOM_WIDTH} --pair-weight ${PAIR_W} \
    --multipion ${HT}/htof6_${cur}.root ${ZERO_LAYERS:+--zero-layers ${ZERO_LAYERS}} \
    > ${CF}/${next}.log 2>&1
  log "fit ${cur} -> ${next}: $(grep -E 'multipion-filled|tracks=' ${CF}/${next}.log | tr '\n' ' ')"
  start_next "${next}"
  log "submitted beam-through + multi-pion Helix (and K-p checks) for ${next}"
  log "QA up to ${cur}:"
  qa "${NAME}_upto_${cur}" "${labels[@]}"
  labels+=("${next}")
  cur="${next}"
done
log "last map ${cur}: wait Helix for final QA"
wait_bt "${cur}"
wait_htof "${cur}"
qa "${NAME}_upto_${cur}" "${labels[@]}"
log "done"
