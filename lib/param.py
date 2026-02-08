from pathlib import Path
import shutil
import sys

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

    def write(self, target_path: Path = None):
        """Write records back to file."""
        if target_path is None:
            target_path = self.path
        
        with open(target_path, 'w') as f:
            for is_comment, parts in self.records:
                if is_comment:
                    # parts[0] is the raw line text for comments
                    f.write(parts[0] + "\n")
                else:
                    # Use tabs for data lines as expected by some scripts
                    f.write("\t".join(parts) + "\n")

    @staticmethod
    def create_from_template(template_path: Path, target_path: Path):
        if not target_path.exists():
            target_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(template_path, target_path)
        return ParamFile(target_path)
