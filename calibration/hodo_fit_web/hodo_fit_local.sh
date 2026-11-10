#!/bin/bash
# Own-PC: ~/hodofit.sh [run] [--suffix Pi]
port="${HODOFIT_PORT:-8767}"
ssh -x -t -L "${port}:localhost:${port}" kek \
  "/home/had/sryuta/JPARC2025E72/calibration/hodo_fit_web/hodo_fit_remote.sh $* --port ${port}" |
while IFS= read -r line; do
  printf '%s\n' "$line"
  case "$line" in *"open http"*) explorer.exe "http://localhost:${port}" >/dev/null 2>&1 ;; esac
done
