"""Atomic JSON publication and reentrant process/thread locks for shared state."""
import fcntl
import json
import os
import tempfile
import threading
from contextlib import contextmanager
from pathlib import Path

_locks = {}
_locks_guard = threading.Lock()
_held = threading.local()


#_____________________________________________________________________________
@contextmanager
def locked_path(path):
    """Serialize a state transaction across threads and hosts using a sidecar."""
    path = Path(path)
    key = str(path.resolve())
    with _locks_guard:
        thread_lock = _locks.setdefault(key, threading.RLock())
    with thread_lock:
        held = getattr(_held, "paths", None)
        if held is None:
            held = _held.paths = set()
        if key in held:
            yield
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        with Path(str(path) + ".lock").open("a") as lock_file:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
            held.add(key)
            try:
                yield
            finally:
                held.remove(key)
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


#_____________________________________________________________________________
def write_json(path, payload):
    """Publish complete JSON with rename on the same filesystem."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    staged = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", delete=False,
        ) as stream:
            staged = Path(stream.name)
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(staged, path)
    finally:
        if staged is not None and staged.exists():
            staged.unlink()
