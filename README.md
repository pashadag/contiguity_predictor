# Contiguity predictor

`contiguity_predictor.py` predicts mean assembly contig length, selects the best
candidate k, and computes the read-arrival rate needed to reach a target fraction
of the limiting mean length. It implements Formula 3, Equation 1, and Lemma 2 in
the manuscript *Theoretical prediction of assembly contiguity*.

Requires Python 3 and uses only the standard library. No packages need installing.
Run the examples below from the directory containing `contiguity_predictor.py`
and this README (the GitHub repository root).

## Parameters

| Parameter | Meaning | Allowed values |
| --- | --- | --- |
| `L` / `--L` | Read length in bases | Integer greater than 1; CLI default: 150 |
| `k` / `--k` | k-mer length, or a list of candidate lengths | Distinct integers from 1 through `L`, inclusive |
| `rho` / `--rho` | Read-arrival rate per base | Finite, nonnegative real number |
| `mu` / `--mu` | k-reftig break rate per base, or corresponding rates for candidate k values | Finite, positive real numbers |
| `sigma` / `--sigma` | Target fraction of the limiting mean length | Real number strictly between 0 and 1 |

For each k, `mu_k` is the reciprocal of the mean reftig length at that k.
For example, a mean reftig length of 1000 bases corresponds to `mu_k = 1e-03`.
The adjusted read length used by the formula is `ell = L - k + 1`.

Convert depth of coverage to a read-arrival rate using `rho = Cov / L`.
For example, coverage 30 with `L = 150` gives `rho = 0.20`.
Convert a returned saturation rate back to depth of coverage using
`Cov_star = rho_star * L`.

## Command-line usage

Place `zval`, `bestk`, or `saturation` immediately after the script name.
Each command prints one result to standard output. Invalid arguments produce an
error on standard error and a nonzero exit status.

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

`--k`, `--rho`, and `--mu` are required. `--L` is optional.
At `rho = 0`, the function returns the continuous limit
`(1 - exp(-mu * L)) / mu`.

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

`--rho`, `--k`, and `--mu` are required. Supply space-separated lists for `--k`
and `--mu`; the lists must be nonempty and equally long. Entries are paired by
position: the example pairs k values 21, 31, and 51 with their respective rates.
Candidates may appear in any order. Ties return the smallest k.

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

`--k`, `--mu`, and `--sigma` are required. The candidate-list rules are the same
as for `bestk`.

The target mean length is `sigma / min(mu_values)`. For the full model sequence
`mu_1 >= ... >= mu_L`, this is `sigma / mu_L`. When supplying only a subset of k
values, the target is relative to the best limiting mean within that subset.

The result is the smallest rate at which at least one candidate reaches the
target. Candidates that cannot reach the target at any finite rate are skipped.
The function returns zero if the zero-rate limit already meets the target.
The implementation evaluates an equivalent logarithmic form of Lemma 2's
Lambert-W expression, without requiring SciPy.

### Help and output precision

```bash
python3 contiguity_predictor.py --help
python3 contiguity_predictor.py zval --help
python3 contiguity_predictor.py bestk --help
python3 contiguity_predictor.py saturation --help
```

`bestk` prints an integer. The other commands print results with two decimal
places, using scientific notation for nonzero magnitudes below `0.01`.
Use the Python functions when retaining full precision for subsequent
calculations, such as converting the saturation rate to coverage.

## Python usage

Import the functions from the directory containing the script:

```python
from contiguity_predictor import zval, bestk, saturation

L = 150
rho = 30 / L
k_values = [21, 31, 51]
mu_values = [3e-03, 2e-03, 1e-03]

mean_length = zval(L, 31, rho, 2e-03)
k_star = bestk(L, rho, k_values, mu_values)
rho_star = saturation(L, k_values, mu_values, 0.95)
coverage_star = rho_star * L

print(f"Predicted mean length: {mean_length:.2f} bases")
print(f"Best k: {k_star}")
print(f"Saturation coverage: {coverage_star:.2f}")
```

The public signatures are:

```python
zval(L, k, rho, mu)                         # Returns a float: mean length in bases.
bestk(L, rho, k_values, mu_values)           # Returns an int: the selected k.
saturation(L, k_values, mu_values, sigma)    # Returns a float: saturation rate rho*.
```

All arguments, including `L`, are required in the Python API. Candidate sequences
are paired by position, as on the command line. The functions return unrounded
values and raise `ValueError` for invalid numeric values or candidate lists.
