from pathlib import Path

# Base paths
# Please modify these paths according to your environment
ANALYZER_DIR = Path("/home/had/sryuta/analyzer/JPARC2025E72")
PARAM_DIR = Path("/home/had/sryuta/analyzer/JPARC2025E72/param")
DATA_DIR = Path("/group/had/sks/Users/sryuta/JPARC2025E72")
OUTPUT_DIR = Path("/home/had/sryuta/JPARC2025E72/results")

# Heavy decoded files (DATA_DIR side): FINAL stages and D5 in DECODE_DIR, others in SCRATCH_DIR
SCRATCH_DIR = DATA_DIR / "scratch_root"
DECODE_DIR = DATA_DIR / "root"

# Backups of params before they are overwritten, and scratch work (numbers, throw-away copies)
BACKUP_DIR = Path("/group/had/sks/Users/sryuta/backup")
WORK_DIR = Path("/group/had/sks/Users/sryuta/tmp")
