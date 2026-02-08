#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Create directories
mkdir -p $SCRIPT_DIR/common/include
mkdir -p $SCRIPT_DIR/lib

# Copy example paths for Python if target doesn't exist
TARGET_FILE="$SCRIPT_DIR/lib/paths.py"
SOURCE_FILE="$SCRIPT_DIR/lib/paths_example.py"

if [ ! -e "$TARGET_FILE" ]; then
    echo "Copying $SOURCE_FILE to $TARGET_FILE"
    cp "$SOURCE_FILE" "$TARGET_FILE"
fi

# Copy example paths for C++ if target doesn't exist
TARGET_FILE_H="$SCRIPT_DIR/common/include/paths.h"
SOURCE_FILE_H="$SCRIPT_DIR/common/include/paths_example.h"

if [ ! -e "$TARGET_FILE_H" ]; then
    echo "Copying $SOURCE_FILE_H to $TARGET_FILE_H"
    cp "$SOURCE_FILE_H" "$TARGET_FILE_H"
fi

# Create output directories
mkdir -p results/root
mkdir -p results/img
mkdir -p data
