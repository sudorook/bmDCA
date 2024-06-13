#!/bin/bash
set -eu

source ../globals

check_dependencies \
  "${PROCESS_MSA_BIN}" \
  "${EXAMPLE_MSA}"

rm -f num_msa.txt weights.txt

"${PROCESS_MSA_BIN}"  -i "${EXAMPLE_MSA}" -g 0.2 -G 0.2 -t 0.8 -o num_msa.txt -O weights.txt

if ! cmp -s num_msa.txt correct_numerical.txt; then
  exit 1
fi

if ! cmp -s weights.txt correct_weights.txt; then
  exit 1
fi
