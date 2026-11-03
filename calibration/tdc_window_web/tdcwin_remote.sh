#!/bin/bash
# KEKCC side of tdcwin_local.sh: start tdc_window_web.py with the node-local envfs python
# (fast imports); fall back to the real conda env. Usage: tdcwin_remote.sh [run] [--trial] [--port N]
"$HOME/local/bin/envfs.sh" mount work >/dev/null 2>&1
py="${XDG_RUNTIME_DIR:+$XDG_RUNTIME_DIR/envfs/work/bin/python3}"
[ -x "$py" ] || py="$HOME/.conda/envs/work/bin/python3"
echo "node $(hostname -s), python $py"
exec "$py" -u "$(dirname "$(readlink -f "$0")")/tdc_window_web.py" "$@"
