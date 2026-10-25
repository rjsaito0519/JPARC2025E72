#!/bin/bash

# List of run numbers to process
RUN_NUMBERS=(
    2140
    2293
    2610
    2775
    2799
    2856
    2884
    2977
    2993
    3005
    3030
)

# Base directory of the analyzer
BASE_DIR="/home/had/sryuta/analyzer/JPARC2025E72/myanalysis"
BIN_EXE="${BASE_DIR}/bin/KVC_analysis"

# Check if the executable exists
if [ ! -f "$BIN_EXE" ]; then
    echo "Error: $BIN_EXE not found. Please run ./build.sh first."
    exit 1
fi

echo "Starting KVC analysis for ${#RUN_NUMBERS[@]} runs..."

for run in "${RUN_NUMBERS[@]}"; do
    echo "--------------------------------------------------"
    echo "Processing Run: $run"
    echo "--------------------------------------------------"
    
    # Run the analysis
    "$BIN_EXE" "$run"
    
    # Check the result
    if [ $? -eq 0 ]; then
        echo "Successfully finished Run: $run"
    else
        echo "Error: Run $run failed."
    fi
done

echo "--------------------------------------------------"
echo "All runs processed."
