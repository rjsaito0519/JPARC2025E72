#!/bin/bash

# ANSI color codes
RED="\033[1;31m"
GREEN="\033[1;32m"
RESET="\033[0m"

# クリーンアップ対象ディレクトリ
BUILD_DIR="./.build"

# クリーン処理を実行
echo "Cleaning build directory..."
rm -rf "${BUILD_DIR}"
if [ $? -eq 0 ]; then
    echo -e "${GREEN}[success]${RESET} Build directory cleaned."
else
    echo -e "${RED}[error]${RESET} Failed to clean build directory."
    exit 1
fi
