#!/bin/bash
set -eu

for dir in *; do
  if [ -d "${dir}" ]; then
    pushd "${dir}" > /dev/null
    if [ -f "run_test.sh" ]; then
      echo -n "${dir}:" >&2
      ./run_test.sh > /dev/null 2>&1
      echo " good."
    fi
    popd > /dev/null
  fi
done
