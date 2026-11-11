import datetime
import fcntl
import shutil
import sys
import time
from contextlib import contextmanager
from pathlib import Path

from . import config


#_____________________________________________________________________________
@contextmanager
def exclusive_param_lock(path, *, timeout_s=600.0):
    """Exclusive flock for read-modify-write of a shared param file.

    Sidecar ``<path>.lock`` so writers can lock before the target exists.
    BcIn/BcOut lanes update the same DCTDC/DCGEO; without this lock one
    chamber's row changes are lost when both do load→edit→save.
    """
    path = Path(path)
    lock_path = Path(str(path) + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_file = lock_path.open("a+", encoding="utf-8")
    deadline = None
    if timeout_s is not None and float(timeout_s) > 0:
        deadline = time.monotonic() + float(timeout_s)
    try:
        while True:
            try:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if deadline is not None and time.monotonic() >= deadline:
                    raise TimeoutError(
                        f"param lock timeout ({timeout_s}s): {lock_path}"
                    )
                time.sleep(0.05)
        yield lock_path
    finally:
        try:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
        finally:
            lock_file.close()


def backup_file(path):
    """Copy path to BACKUP_DIR/param_history/<date>/<name>.<time> before it is overwritten; the copy, or None.
    (Not BACKUP_DIR/param: that directory is an rsync --delete mirror of param/ and loses anything else.)"""
    path = Path(path)
    if not path.exists():
        return None
    now = datetime.datetime.now()
    dst = config.BACKUP_DIR / "param_history" / f"{now:%Y%m%d}" / f"{path.name}.{now:%H%M%S}.{now.microsecond // 1000:03d}"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, dst)
    print(f"backup: {dst}")
    return dst


class ParamFile:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.records = [] # List of (is_comment, parts_list)
        if self.path.exists():
            self._read()

    def _read(self):
        self.records = []
        with open(self.path, 'r') as f:
            for line in f:
                s_list = line.split()
                if not s_list or s_list[0].startswith('#'):
                    # Keep raw line text for comments to preserve layout
                    self.records.append((True, [line.strip()]))
                else:
                    self.records.append((False, s_list))

    def get_data_dict(self, key_len=1):
        """
        Get current data as a dictionary.
        key_len: Number of columns to join as the key (default 1).
        """
        data = {}
        for is_comment, parts in self.records:
            if not is_comment and len(parts) >= key_len:
                key = "-".join(parts[:key_len])
                data[key] = parts[key_len:]
        return data

    def key_issues(self, new_data, key_len=1):
        """Compare new_data keys to this file; return (missing, duplicates).

        missing: keys in new_data with no matching data line
        duplicates: keys in new_data that appear more than once in the file
        Does not change records; ParamFile.update keeps its silent skip behavior.
        """
        counts = {}
        for is_comment, parts in self.records:
            if not is_comment and len(parts) >= key_len:
                key = "-".join(parts[:key_len])
                counts[key] = counts.get(key, 0) + 1
        missing = sorted(k for k in new_data if counts.get(k, 0) == 0)
        duplicates = sorted(k for k in new_data if counts.get(k, 0) > 1)
        return missing, duplicates

    def update(self, new_data, key_len=1, start_col=None):
        """
        Surgically update existing data with new values.
        new_data: {key: [val1, val2, ...]}
        key_len: Number of columns that make up the key.
        start_col: Index where the values in new_data[key] should be placed.
                   Defaults to key_len (immediately after the key).
        """
        if start_col is None:
            start_col = key_len
            
        updated_count = 0
        for i, (is_comment, parts) in enumerate(self.records):
            if not is_comment and len(parts) >= key_len:
                key = "-".join(parts[:key_len])
                if key in new_data:
                    old_parts = list(parts)
                    for j, v in enumerate(new_data[key]):
                        target_idx = start_col + j
                        if target_idx < len(old_parts):
                            old_parts[target_idx] = str(v)
                        else:
                            # If new data has more values than existing columns, extend.
                            old_parts.append(str(v))
                    self.records[i] = (False, old_parts)
                    updated_count += 1
        return updated_count

    def serialized_text(self):
        """Same byte layout as write(), without touching the file (for content hashing)."""
        lines = []
        for is_comment, parts in self.records:
            if is_comment:
                lines.append(parts[0] + "\n")
            else:
                lines.append("\t".join(parts) + "\n")
        return "".join(lines)

    def write(self, target_path: Path = None):
        """Write records back to file."""
        if target_path is None:
            target_path = self.path
        backup_file(target_path)
        with open(target_path, 'w') as f:
            f.write(self.serialized_text())

    @staticmethod
    def create_from_template(template_path: Path, target_path: Path):
        if not target_path.exists():
            target_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(template_path, target_path)
        return ParamFile(target_path)
