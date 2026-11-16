#!/bin/bash
# WSL launcher for Hodoscope calibration review (ADC and T0): ~/calib.sh RUN [--kaon].
set -euo pipefail
if [ "$#" -lt 1 ] || [[ ! "$1" =~ ^[0-9]{1,5}$ ]]; then
  echo 'Usage: calib.sh RUN [--kaon]' >&2
  exit 1
fi
run="$1"
shift
particle_arg=""
if [ "$#" -gt 0 ]; then
  if [ "$#" -ne 1 ] || [ "$1" != "--kaon" ]; then
    echo 'Usage: calib.sh RUN [--kaon]' >&2
    exit 1
  fi
  particle_arg="--kaon"
fi
port="${CALIB_PORT:-${MIP_PORT:-8767}}"
if [[ ! "$port" =~ ^[0-9]{4,5}$ ]] || [ "$port" -lt 1024 ] || [ "$port" -gt 65535 ]; then
  echo 'CALIB_PORT must be a port number (1024–65535)' >&2
  exit 1
fi
ssh -x -t -o ExitOnForwardFailure=yes -L "${port}:localhost:${port}" kek \
  "/home/had/sryuta/JPARC2025E72/calibration/mip_web/mip_remote.sh ${run} ${particle_arg} --port ${port}" |
while IFS= read -r line; do
  printf '%s\n' "$line"
  case "$line" in
    *"open http"*)
      if command -v explorer.exe >/dev/null 2>&1; then
        explorer.exe "http://localhost:${port}/?run=${run}" >/dev/null 2>&1 || true
      else
        echo "Browser: http://localhost:${port}/?run=${run}"
      fi ;;
  esac
done
