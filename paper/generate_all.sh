#!/usr/bin/env bash
# Regenerate Table S1, Figure 5 and Figure 6 of the paper with default parameters.
set -euo pipefail

cd "$(dirname "$0")"

echo "Generating Table S1 (suppl_table.csv)..."
python3 generate_s1_table.py

echo "Generating Figure 5 (figure_5.pdf)..."
python3 figure_5.py

echo "Generating Figure 6 (figure_6.pdf)..."
python3 figure_6.py

echo "Done. Outputs are in $(pwd)"
