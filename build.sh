#!/bin/bash

# ANSI color codes
RED="\033[1;31m"
GREEN="\033[1;32m"
RESET="\033[0m"

# ビルドディレクトリの指定
BUILD_DIR="./.build"

# 必要なディレクトリを作成
echo "Creating build directory..."
mkdir -p "${BUILD_DIR}"
if [ $? -ne 0 ]; then
    echo -e "${RED}[error]${RESET} Failed to create build directory."
    exit 1
fi
echo -e "${GREEN}[success]${RESET} Build directory created."

# CMakeを実行
echo "Running CMake configuration..."
cd "${BUILD_DIR}"
cmake .. 2>&1 | tee cmake_config.log
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo -e "${RED}[error]${RESET} CMake configuration failed! Check ${BUILD_DIR}/cmake_config.log for details."
    exit 1
fi
echo -e "${GREEN}[success]${RESET} CMake configuration completed."

# Makeを実行
echo "Building the project with half the available CPU cores..."
make -j$(($(nproc)/2)) 2>&1 | tee build.log
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo -e "${RED}[error]${RESET} Build failed! Check ${BUILD_DIR}/build.log for details."
    exit 1
fi
echo -e "${GREEN}[success]${RESET} Build successful!"
cd ..