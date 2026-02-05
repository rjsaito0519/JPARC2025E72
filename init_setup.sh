#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir $SCRIPT_DIR/results
mkdir $SCRIPT_DIR/results/root
mkdir $SCRIPT_DIR/results/img

mkdir $SCRIPT_DIR/data

TARGET_FILE="$SCRIPT_DIR/include/paths.h"
SOURCE_FILE="$SCRIPT_DIR/include/paths_example.h"

if [ ! -e "$TARGET_FILE" ]; then
  cp "$SOURCE_FILE" "$TARGET_FILE"
fi

TARGET_FILE="$SCRIPT_DIR/scripts/conf.py"
SOURCE_FILE="$SCRIPT_DIR/scripts/conf_example.py"

if [ ! -e "$TARGET_FILE" ]; then
  cp "$SOURCE_FILE" "$TARGET_FILE"
fi
