#!/bin/bash
set -eu

source ../globals

check_dependencies

rm -rf results_continuous results_2_steps
rm -f num_msa.txt weights.txt

"${PROCESS_MSA_BIN}"  -i "${EXAMPLE_MSA}" -g 0.2 -G 0.2 -t 0.8 -o num_msa.txt -O weights.txt

"${BMDCA_BIN}" -n num_msa.txt -w weights.txt -c <(sed -e "s/^step_max=.*/step_max=10/g" -e "s/^save_period=.*/save_period=10/g" "${EXAMPLE_BMDCA_CONF}")  -d results_continuous

"${BMDCA_BIN}" -n num_msa.txt -w weights.txt -c <(sed -e "s/^step_max=.*/step_max=5/g" -e "s/^save_period=.*/save_period=5/g" "${EXAMPLE_BMDCA_CONF}")  -d results_2_steps
"${BMDCA_BIN}" -n num_msa.txt -w weights.txt -c <(sed -e "s/^step_max=.*/step_max=10/g" -e "s/^save_period=.*/save_period=5/g" "${EXAMPLE_BMDCA_CONF}")  -d results_2_steps

for FILE in results_continuous/*{.bin,.txt}; do
  if ! cmp -s "${FILE}" "${FILE//continuous/2_steps}"; then
    echo "${FILE@Q} mismatch."
    exit 1
  fi
done
