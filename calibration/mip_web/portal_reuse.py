"""Keep an SSH connection open when this node already runs a calibration portal."""
import argparse
import json
import time
import urllib.request
parser = argparse.ArgumentParser()
parser.add_argument("run", type=int, nargs="?")
parser.add_argument("--kaon", action="store_true")
parser.add_argument("--port", type=int, default=8767)
args = parser.parse_args()
try:
    with urllib.request.urlopen(f"http://127.0.0.1:{args.port}/portal-health", timeout=2) as response:
        health = json.load(response)
    if health != {"service": "e72_calib_portal", "version": 1}:
        raise ValueError("different service")
except (OSError, ValueError):
    raise SystemExit(1)
query = f"?run={args.run}" if args.run is not None else ""
if args.kaon:
    query += ("&" if query else "?") + "particle=K"
print(f"open http://localhost:{args.port}/{query}  (existing calibration portal; Ctrl-C closes connection)", flush=True)
try:
    while True:
        time.sleep(30)
except KeyboardInterrupt:
    pass
