#!/bin/bash
set -e
"$HOME/local/bin/envfs.sh" mount work >/dev/null 2>&1 || true
python_work="${XDG_RUNTIME_DIR:+$XDG_RUNTIME_DIR/envfs/work/bin/python3}"
[ -x "$python_work" ] || python_work="$HOME/.conda/envs/work/bin/python3"
echo "node $(hostname -s), conda work: $python_work"
set +e
"$python_work" -u "$(dirname "$(readlink -f "$0")")/portal_reuse.py" "$@"
reuse_status=$?
set -e
if [ "$reuse_status" -eq 0 ]; then exit 0; fi
if [ "$reuse_status" -ne 1 ]; then exit "$reuse_status"; fi
exec "$python_work" -u "$(dirname "$(readlink -f "$0")")/mip_web.py" "$@"
