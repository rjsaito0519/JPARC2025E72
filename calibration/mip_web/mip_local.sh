#!/bin/bash
# Compatibility entry point for older installations.
exec "$(dirname "$(readlink -f "$0")")/calib_local.sh" "$@"
