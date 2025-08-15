#!/bin/bash

experiment=$1
i=$2

if [ -z "$experiment" ] || [ -z "$i" ]; then
  echo "Usage: $0 <experiment> <job_index>"
  echo "Example: $0 testSEI3R_02 3"
  exit 1
fi

# Get absolute path to the directory containing this script
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Define the output directory
output_dir="$BASE_DIR/haics/simulation/output/${experiment}_$i"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Change to the correct working directory (where ./bin/haics is located)
cd "$BASE_DIR/haics/simulation" || exit 1

# Run with time profiling
/usr/bin/time -f "Elapsed: %E | CPU: %P | Mem: %M KB" \
  ./bin/haics "${experiment}_$i" \
  > "${output_dir}/log.out" \
  2> "${output_dir}/log.err"

echo "Finished job $i" >> "${output_dir}/log.out"
