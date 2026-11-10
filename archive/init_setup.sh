#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$SCRIPT_DIR/common/include"
mkdir -p "$SCRIPT_DIR/common/lib"

TARGET_FILE="$SCRIPT_DIR/common/lib/paths.py"
SOURCE_FILE="$SCRIPT_DIR/common/lib/paths_example.py"

if [ ! -e "$TARGET_FILE" ]; then
    if [ ! -f "$SOURCE_FILE" ]; then
        echo "Error: missing $SOURCE_FILE" >&2
        exit 1
    fi
    echo "Copying $SOURCE_FILE to $TARGET_FILE"
    cp "$SOURCE_FILE" "$TARGET_FILE"
fi

TARGET_FILE_H="$SCRIPT_DIR/common/include/paths.h"
SOURCE_FILE_H="$SCRIPT_DIR/common/include/paths_example.h"

if [ ! -e "$TARGET_FILE_H" ]; then
    if [ ! -f "$SOURCE_FILE_H" ]; then
        echo "Error: missing $SOURCE_FILE_H" >&2
        exit 1
    fi
    echo "Copying $SOURCE_FILE_H to $TARGET_FILE_H"
    cp "$SOURCE_FILE_H" "$TARGET_FILE_H"
fi

mkdir -p results/root
mkdir -p results/img
mkdir -p data
