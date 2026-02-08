import time
import sys
from pathlib import Path

start = time.time()
print(f"Start: {time.time() - start:.4f}s")

import argparse
print(f"argparse: {time.time() - start:.4f}s")
import numpy as np
print(f"numpy: {time.time() - start:.4f}s")
import matplotlib.pyplot as plt
print(f"matplotlib: {time.time() - start:.4f}s")
import uproot
print(f"uproot: {time.time() - start:.4f}s")

project_root = Path("/home/had/sryuta/analyzer/JPARC2025E72/myanalysis")
sys.path.append(str(project_root))
from lib import config
print(f"config: {time.time() - start:.4f}s")
