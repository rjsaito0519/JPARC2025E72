#!/bin/bash
# Own-PC (WSL) side. Copy to ~/tdcwin.sh on the PC, then: ~/tdcwin.sh [run] [--trial]
# Opens an ssh tunnel to KEKCC (ssh host alias "kek"), starts the server there and opens the
# page in the Windows browser once it is ready. Ctrl-C stops both.
port="${TDCWIN_PORT:-8765}"
ssh -x -t -L "${port}:localhost:${port}" kek \
  "/home/had/sryuta/JPARC2025E72/calibration/tdc_window_web/tdcwin_remote.sh $* --port ${port}" |
while IFS= read -r line; do
  printf '%s\n' "$line"
  case "$line" in *"open http"*) explorer.exe "http://localhost:${port}" >/dev/null 2>&1 ;; esac
done
