# Reproducing the paper's results

This folder regenerates Table S1, Figure 5 and Figure 6 of our paper *Theoretical prediction of assembly contiguity* from the data files in this directory.

## Requirements

Python 3 with SciPy, NumPy and Matplotlib:

```bash
python3 -m pip install scipy numpy matplotlib
```

## Usage

Run the following script from `paper` directory.

```bash
bash generate_all.sh
```

It runs `generate_s1_table.py`, `figure_5.py` and `figure_6.py` with their default parameters.

## Expected output

Following files are generated inside the `paper` folder:

| File | Contents |
|------|----------|
| `suppl_table_s1.csv` | Table S1 |
| `figure_5.pdf` | Figure 5 |
| `figure_6.pdf` | Figure 6 |
