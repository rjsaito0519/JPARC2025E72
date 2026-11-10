#!/bin/bash
"$HOME/local/bin/envfs.sh" mount work >/dev/null 2>&1
py="${XDG_RUNTIME_DIR:+$XDG_RUNTIME_DIR/envfs/work/bin/python3}"
[ -x "$py" ] || py="$HOME/.conda/envs/work/bin/python3"
echo "node $(hostname -s), python $py"
exec "$py" -u "$(dirname "$(readlink -f "$0")")/hodo_fit_web.py" "$@"
