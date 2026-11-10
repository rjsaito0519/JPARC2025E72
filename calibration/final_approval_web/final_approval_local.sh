#!/bin/bash
# Own-PC (WSL) side. Copy to ~/finalqa.sh then: ~/finalqa.sh [run]
port="${FINALQA_PORT:-8766}"
ssh -x -t -L "${port}:localhost:${port}" kek \
  "/home/had/sryuta/JPARC2025E72/calibration/final_approval_web/final_approval_remote.sh $* --port ${port}" |
while IFS= read -r line; do
  printf '%s\n' "$line"
  case "$line" in *"open http"*) explorer.exe "http://localhost:${port}" >/dev/null 2>&1 ;; esac
done
