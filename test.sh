#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_command="${PYTHON:-python3}"

# P. brunnea (GCF_037157775.1): rounded reftig-break rates for k = 20, 30, ..., 140.
k_values=(20 30 40 50 60 70 80 90 100 110 120 130 140)
mu_values=(
    3.93e-03 4.94e-04 2.06e-04 1.14e-04 7.72e-05 5.81e-05 4.38e-05
    3.63e-05 2.95e-05 2.48e-05 2.13e-05 1.86e-05 1.65e-05
)

check() {
    local name="$1"
    local expected="$2"
    shift 2
    local actual
    actual=$("$python_command" "$script_dir/contiguity_predictor.py" "$@")
    if [[ "$actual" != "$expected" ]]; then
        printf 'FAIL %s: expected %s, got %s\n' "$name" "$expected" "$actual" >&2
        exit 1
    fi
    printf 'PASS %s: %s\n' "$name" "$actual"
}

check zval 26809.13 zval \
    --L 150 --k 90 --rho 0.20 --mu "${mu_values[7]}"

check bestk 100 bestk \
    --L 150 --rho 0.20 --k "${k_values[@]}" --mu "${mu_values[@]}"

check saturation 1.29 saturation \
    --L 150 --k "${k_values[@]}" --mu "${mu_values[@]}" --sigma 0.95

printf 'All installation checks passed.\n'
