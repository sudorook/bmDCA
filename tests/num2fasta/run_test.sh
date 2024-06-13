#!/bin/bash
set -eu

source ../globals

check_dependencies \
  "${NUMERIC2FASTA_BIN}"

rm -f tmp.fasta
"${NUMERIC2FASTA_BIN}" -n source.txt -o tmp.fasta

# should pass
if ! cmp -s tmp.fasta target.fasta; then
  exit 1
fi

"${NUMERIC2FASTA_BIN}" -n <(head -n 100 source.txt) -o tmp.fasta

# should fail
if cmp -s tmp.fasta target.fasta; then
  exit 1
fi
