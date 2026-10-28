# Contiguity predictor

`contiguity_predictor.py` predicts mean assembly contig length, selects the best
candidate k, and computes the read-arrival rate needed to reach a target fraction
of the limiting mean length. It implements Formula 3, Equation 1, and Lemma 2 in
the manuscript *Theoretical prediction of assembly contiguity*.

Requires Python 3 and SciPy. Install SciPy with:

```bash
python3 -m pip install scipy
```

Run the examples below from the directory containing `contiguity_predictor.py`
and this README (the GitHub repository root).

## Verify the installation

Run:

```bash
bash test.sh
```

This checks all three commands using hardcoded P. brunnea rates.
Each command should report `PASS`.
To select a different Python interpreter, use `PYTHON=/path/to/python bash test.sh`.

## Parameters

| Parameter | Meaning | Allowed values |
| --- | --- | --- |
| `L` / `--L` | Read length in bases | Integer greater than 1; required |
| `k` / `--k` | k-mer length, or a list of candidate lengths | Distinct integers from 1 through `L`, inclusive |
| `rho` / `--rho` | Read-arrival rate per base | Real number greater than 0 |
| `mu` / `--mu` | k-reftig break rate per base, or corresponding rates for candidate k values | Real numbers greater than 0; non-increasing as k increases for `bestk` and `saturation` |
| `sigma` / `--sigma` | Target fraction of the limiting mean length | Real number strictly between 0 and 1 |

## Command-line usage

### Predict mean contig length: `zval`

Evaluate Formula 3 for one k and its break rate:

```bash
python3 contiguity_predictor.py zval \
    --L 150 --k 31 --rho 0.20 --mu 1e-03
```

Output, in bases:

```text
1000.00
```

`--L`, `--k`, `--rho`, and `--mu` are required.

### Select the best k: `bestk`

Evaluate Equation 1 over the supplied candidate k values:

```bash
python3 contiguity_predictor.py bestk \
    --L 150 --rho 0.20 --k 21 31 51 --mu 3e-03 2e-03 1e-03
```

Output:

```text
51
```

`--L`, `--rho`, `--k`, and `--mu` are required. Supply space-separated lists for `--k`
and `--mu`; the lists must be nonempty and equally long. Entries are paired by
position: the example pairs k values 21, 31, and 51 with their respective rates.
Candidates may appear in any order. Their corresponding `mu_k` values must be
non-increasing as k increases; equal values are allowed. Ties return the smallest k.

### Find the saturation rate: `saturation`

Use Lemma 2 to compute the minimum read-arrival rate over the supplied candidates:

```bash
python3 contiguity_predictor.py saturation \
    --L 150 --k 21 31 51 --mu 3e-03 2e-03 1e-03 --sigma 0.95
```

Output, as a read-arrival rate:

```text
0.07
```

`--L`, `--k`, `--mu`, and `--sigma` are required. The candidate-list rules are the same
as for `bestk`.

The target mean length is `sigma / min(mu_values)`. For the full model sequence
`mu_1 >= ... >= mu_L`, this is `sigma / mu_L`. When supplying only a subset of k
values, the target is relative to the best limiting mean within that subset.

The result is the smallest rate at which at least one candidate reaches the
target.

## Python usage

Import the functions from the directory containing the script:

```python
from contiguity_predictor import zval, bestk, saturation

L = 150
rho = 0.20
k_values = [21, 31, 51]
mu_values = [3e-03, 2e-03, 1e-03]

mean_length = zval(L, 31, rho, 2e-03)
k_star = bestk(L, rho, k_values, mu_values)
rho_star = saturation(L, k_values, mu_values, 0.95)

print(f"Predicted mean length: {mean_length:.2f} bases")
print(f"Best k: {k_star}")
print(f"Saturation rate: {rho_star:.2f}")
```
