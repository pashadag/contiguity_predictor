# Reproducing the paper's results

This folder regenerates Table S1, Figure 5 and Figure 6 of our paper *Theoretical prediction of assembly contiguity* from the data files in this directory.

## Requirements

Python 3 with SciPy, NumPy and Matplotlib:

```bash
python3 -m pip install scipy numpy matplotlib
```

## Usage

```bash
bash generate_all.sh
```

The script can be run from any directory; it runs `generate_s1_table.py`, `figure_5.py` and `figure_6.py` with their default parameters.

## Expected output

The script prints a line for each file as it is generated:

```
Generating Table S1 (suppl_table.csv)...
wrote 84 rows to .../paper/suppl_table.csv
Generating Figure 5 (figure_5.pdf)...
gamma = 1: no crossing
gamma = 5: M2 coverage = 3.14, M1 coverage = 0.63, Z = 370
gamma = 5: M2 coverage = 71.71, M1 coverage = 14.34, Z = 33,704
gamma = 10: M2 coverage = 2.82, M1 coverage = 0.28, Z = 326
gamma = 10: M2 coverage = 157.37, M1 coverage = 15.74, Z = 40,250
Plot saved to .../paper/figure_5.pdf
Generating Figure 6 (figure_6.pdf)...
Saved figure to .../paper/figure_6.pdf
Done. Outputs are in .../paper
```

and writes these files to the `paper` folder:

| File | Contents |
|------|----------|
| `suppl_table.csv` | Table S1 |
| `figure_5.pdf` | Figure 5 |
| `figure_6.pdf` | Figure 6 |
