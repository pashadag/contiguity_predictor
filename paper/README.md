# Reproducing the paper's results

These scripts regenerate Table S1 and Figure 5 of our paper *Theoretical prediction of assembly contiguity* from the data files in this directory. Run them from this folder.

## Requirements: 
Python 3 with SciPy (`generate_table.py`), and NumPy and
Matplotlib (`budget_plot.py`).

```bash
python3 -m pip install scipy numpy matplotlib
```

## Table S1: `generate_table.py`

```bash
python3 generate_table.py [--out suppl_table.csv]
```

The script uses `CSV` files from `paper` folder, `zval` and `bestk` from `contiguity_predictor.py`, with `L = 150`, `rho = coverage / L` and k in {20, 30, ..., 150}.

## Figure 5: `budget_plot.py`

```bash
python3 budget_plot.py [--stats p_brunnea_reftig_stats.txt] \
    [--out two-machine-budget.pdf]
```

The script plots the figure 5 from our paper.
