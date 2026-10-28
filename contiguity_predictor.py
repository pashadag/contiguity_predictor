#!/usr/bin/env python3
"""Predict assembly contiguity (Python standard library only).

Based on the manuscript "Theoretical prediction of assembly contiguity".

Public functions:
    zval(L, k, rho, mu): expected unitig length, Formula 3.
    bestk(L, rho, k_values, mu_values): maximizing candidate k, Equation 1.
    saturation(L, k_values, mu_values, sigma): minimum read-arrival rate,
        using the fixed-k saturation formula in Lemma 2.

rho is the read-arrival rate, not depth of coverage: rho = Cov / L.
For a subset of candidate k values, saturation is relative to the best limiting
mean in that subset, 1 / min(mu_values). For the full, nonincreasing model
sequence, this is 1 / mu_L. Rates and lengths returned by the functions retain
full precision; command-line results use at most two decimal places.

Examples:
    python3 contiguity_predictor.py zval --L 150 --k 31 --rho 0.20 --mu 0.001
    python3 contiguity_predictor.py bestk --L 150 --rho 0.20 --k 21 31 51 --mu 0.003 0.002 0.001
    python3 contiguity_predictor.py saturation --L 150 --k 21 31 51 --mu 0.003 0.002 0.001 --sigma 0.95

The command immediately following the script name selects the function.
Use --help after any command for its arguments; --L is required.
"""

import argparse
import math
from numbers import Integral


__all__ = ["zval", "bestk", "saturation"]


def _read_length(L):
    if isinstance(L, bool) or not isinstance(L, Integral) or L <= 1:
        raise ValueError("L must be an integer greater than 1")
    return int(L)


def _candidate_k(k, L):
    if isinstance(k, bool) or not isinstance(k, Integral) or not 1 <= k <= L:
        raise ValueError("each k must be an integer between 1 and L, inclusive")
    return int(k)


def _finite_number(value, name):
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{name} must be a finite real number") from None
    if isinstance(value, bool) or not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _candidates(L, k_values, mu_values):
    L = _read_length(L)
    ks, mus = list(k_values), list(mu_values)
    if not ks or len(ks) != len(mus):
        raise ValueError("k_values and mu_values must be nonempty and equally long")
    ks = [_candidate_k(k, L) for k in ks]
    if len(set(ks)) != len(ks):
        raise ValueError("candidate k values must be distinct")
    mus = [_finite_number(mu, "mu_k") for mu in mus]
    candidates = list(zip(ks, mus))
    ordered = sorted(candidates)
    if any(left[1] < right[1] for left, right in zip(ordered, ordered[1:])):
        raise ValueError("mu_k values must be non-increasing as k increases")
    return L, candidates


def _softplus(x):
    """log(1 + exp(x)), without exponential overflow."""
    return max(x, 0.0) + math.log1p(math.exp(-abs(x)))


def _zval(L, k, rho, mu):
    ell = L - k + 1
    # Rewrite the fraction in Formula 3 as exp(-mu * (k - 1)) / D,
    # where D = 1 + mu / (rho + mu) * expm1((rho + mu) * ell).
    # This avoids subtracting almost equal numbers when mu is very small.
    exponent = (rho + mu) * ell
    if math.isinf(exponent):
        return 1.0 / mu
    if exponent < 50.0:
        log_d = math.log1p(mu / (rho + mu) * math.expm1(exponent))
    else:
        log_d = _softplus(math.log(mu) - math.log(rho + mu) + exponent)
    return -math.expm1(-mu * (k - 1) - log_d) / mu


def zval(L, k, rho, mu):
    """Return Formula 3's expected unitig length for ell = L - k + 1.

    L must be an integer greater than 1, k an integer in [1, L], rho a
    finite positive read-arrival rate, and mu a finite positive break rate.
    Invalid inputs raise ValueError.
    """
    L = _read_length(L)
    k = _candidate_k(k, L)
    rho = _finite_number(rho, "rho")
    mu = _finite_number(mu, "mu")
    return _zval(L, k, rho, mu)


def bestk(L, rho, k_values, mu_values):
    """Return the candidate k maximizing Formula 3, as in Equation 1.

    Candidate sequences must be nonempty, equally long, and paired by position.
    k values must be distinct integers in [1, L], and break rates positive and
    finite. Break rates must be non-increasing as k increases; equal rates are
    allowed. Candidates may be supplied in any order. Ties return the smallest k.
    rho must be positive and finite. Invalid inputs raise ValueError.
    """
    L, candidates = _candidates(L, k_values, mu_values)
    rho = _finite_number(rho, "rho")
    return max(
        candidates,
        key=lambda pair: (_zval(L, pair[0], rho, pair[1]), -pair[0]),
    )[0]


def _fixed_k_saturation(L, k, mu, reference_mu, sigma):
    ratio = reference_mu / mu
    # A candidate's infinite-rate mean must strictly exceed the target.
    if ratio <= sigma:
        return math.inf
    if ratio * (-math.expm1(-mu * L)) >= sigma:
        return 0.0

    ell = L - k + 1
    # Lemma 2: A = exp(-mu * L) / (1 - sigma * mu / reference_mu),
    # B = (A - exp(-mu * ell)) / mu, with A > 1 and B > 0 here.
    fraction = sigma / ratio
    log_delta = (math.log1p(-fraction) if fraction < 1.0
                 else math.log(ratio - sigma) - math.log(ratio))
    log_a = -mu * L - log_delta
    if log_a <= 0.0:
        return 0.0
    log_d = (math.log(-math.expm1(-mu * ell - log_a))
             - math.log(mu) - math.log(ell))  # d = B / (ell * A)

    # Evaluate the W_{-1} expression through its equivalent equation
    # y = log(A) + log(1 + d * y), y = ell * rho.
    # Its unique positive root selects the -1 branch. Working with logarithms
    # avoids an underflowing Lambert-W argument and requires no SciPy install.
    def residual(y):
        return y - log_a - _softplus(log_d + math.log(y))

    low = log_a
    high = max(1.0, log_a, 2.0 * log_d)
    while residual(high) < 0.0:
        high *= 2.0
    for _ in range(1100):
        midpoint = low + (high - low) / 2.0
        if midpoint == low or midpoint == high or high - low <= 1e-12 * high:
            break
        if residual(midpoint) >= 0.0:
            high = midpoint
        else:
            low = midpoint
    return high / ell


def saturation(L, k_values, mu_values, sigma):
    """Return rho* using Lemma 2, minimized over the supplied candidates.

    The target mean is sigma / min(mu_values); sigma must be finite and strictly
    between 0 and 1. With the full model sequence this is sigma / mu_L. With a
    subset, the reference is the best infinite-rate mean among those candidates.
    Return 0 if the zero-rate limit already meets the target. Candidates whose
    infinite-rate means cannot exceed the target are skipped.

    Candidate validation is the same as bestk. The result is a read-arrival rate;
    multiply by L for depth of coverage. Invalid inputs raise ValueError.
    """
    L, candidates = _candidates(L, k_values, mu_values)
    sigma = _finite_number(sigma, "sigma")
    if sigma >= 1.0:
        raise ValueError("sigma must be strictly between 0 and 1")
    reference_mu = min(mu for _, mu in candidates)
    return min(_fixed_k_saturation(L, k, mu, reference_mu, sigma)
               for k, mu in candidates)


def _format_result(value):
    if value != 0 and abs(value) < 0.01:
        return f"{value:.2e}"
    return f"{value:.2f}"


def _main(argv=None):
    parser = argparse.ArgumentParser(
        description="Predict assembly contiguity; rho = Cov / L.",
        epilog="Example: %(prog)s zval --L 150 --k 31 --rho 0.20 --mu 0.001",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    z_parser = commands.add_parser("zval", help="expected unitig length (Formula 3)")
    k_parser = commands.add_parser("bestk", help="maximizing candidate k (Equation 1)")
    s_parser = commands.add_parser("saturation", help="saturating rho* (Lemma 2)")
    for command in (z_parser, k_parser, s_parser):
        command.add_argument("--L", type=int, required=True, help="read length")
    z_parser.add_argument("--k", type=int, required=True, help="k-mer length")
    z_parser.add_argument("--rho", type=float, required=True, help="positive read-arrival rate")
    z_parser.add_argument("--mu", type=float, required=True, help="positive k-reftig break rate")
    k_parser.add_argument("--rho", type=float, required=True, help="positive read-arrival rate")
    for command in (k_parser, s_parser):
        command.add_argument("--k", type=int, nargs="+", required=True,
                             help="space-separated candidate k values")
        command.add_argument("--mu", type=float, nargs="+", required=True,
                             help="space-separated positive mu_k values, non-increasing as k increases")
    s_parser.add_argument("--sigma", type=float, required=True,
                         help="target fraction, strictly between 0 and 1")
    args = parser.parse_args(argv)
    try:
        if args.command == "zval":
            result = zval(args.L, args.k, args.rho, args.mu)
        elif args.command == "bestk":
            result = bestk(args.L, args.rho, args.k, args.mu)
        else:
            result = saturation(args.L, args.k, args.mu, args.sigma)
    except ValueError as error:
        parser.error(str(error))
    print(result if args.command == "bestk" else _format_result(result))


if __name__ == "__main__":
    _main()
