#!/usr/bin/env python3
"""Generate the accuracy table (Table S1) from the reftig and unitigome stats.

For each genome and coverage c in {5, 9, ..., 25}, with L = 150 and rho = c / L:
    Zk_rho-value error (%): mean absolute percentage error of the Z_{k,rho}-value
        against the Alg. 1 mean unitig length, over k in {20, 30, ..., 150}.
    Z_rho-value (bp): Z_{k*,rho}-value, with k* the predicted best k.
    Mean contig length (bp): maximum Alg. 1 mean unitig length over k.
    Z_rho-value error (%): absolute percentage error of the Z_rho-value
        against the mean contig length.
    predicted k*: k maximizing the Z_{k,rho}-value.
    Alg. 1 opt k: k maximizing the Alg. 1 mean unitig length.

mu_k is 1 / Mean from paper/reftig_stats.csv; Alg. 1 lengths are the Mean
column of paper/unitigome_stats.csv. Values are written at full precision.

Usage:
    python3 generate_table.py [--out generated_table.csv]
"""

import argparse
import csv
import sys
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
ROOT = TEST_DIR.parent
sys.path.insert(0, str(ROOT))

from contiguity_predictor import bestk, zval  

L = 150
K_VALUES = list(range(20, 151, 10))
COVERAGES = list(range(5, 26, 4))

COLUMNS = [
    "Species",
    "Coverage",
    "Zk_rho-value error (%)",
    "Z_rho-value (bp)",
    "Mean contig length (bp)",
    "Z_rho-value error (%)",
    "predicted k*",
    "Alg. 1 opt k",
]

def read_genome_names(path):
    names = {}
    with open(path) as handle:
        for line in handle:
            if line.strip():
                accession, name = line.strip().split(None, 1)
                names[accession] = name
    return names


def read_mu(path):
    """Return {accession: {k: mu_k}} for k in K_VALUES."""
    mu = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            k = int(row["k-mer"])
            if k in K_VALUES:
                mu.setdefault(row["Accession"], {})[k] = 1.0 / float(row["Mean"])
    return mu


def read_unitig_means(path):
    """Return {(accession, coverage): {k: mean unitig length}}."""
    means = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            accession = "_".join(row["Genome"].split("_")[:2])
            k, coverage = int(row["k-mer"]), int(row["Coverage"])
            if k in K_VALUES and coverage in COVERAGES:
                means.setdefault((accession, coverage), {})[k] = float(row["Mean"])
    return means


def percentage_error(predicted, true):
    return abs(predicted - true) / true * 100.0


def table_rows(names, mu, means):
    for accession in names.keys():
        species = names[accession]
        mu_k = [mu[accession][k] for k in K_VALUES]
        for coverage in COVERAGES:
            true = means.get((accession, coverage))
            if true is None or set(true) != set(K_VALUES):
                print(f"skipping {accession} at coverage {coverage}: "
                      "incomplete unitigome stats", file=sys.stderr)
                continue
            rho = coverage / L
            zk_error = sum(
                percentage_error(zval(L, k, rho, mu[accession][k]), true[k])
                for k in K_VALUES
            ) / len(K_VALUES)
            k_star = bestk(L, rho, K_VALUES, mu_k)
            z_rho = zval(L, k_star, rho, mu[accession][k_star])
            # Ties go to the smallest k, as in bestk.
            opt_k = max(K_VALUES, key=lambda k: (true[k], -k))
            mean_contig = true[opt_k]
            yield [
                species,
                coverage,
                zk_error,
                z_rho,
                mean_contig,
                percentage_error(z_rho, mean_contig),
                k_star,
                opt_k,
            ]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=TEST_DIR / "suppl_table_s1.csv")
    args = parser.parse_args()

    names = read_genome_names(TEST_DIR / "genome_names.txt")
    mu = read_mu(ROOT / "paper" / "reftig_stats.csv")
    means = read_unitig_means(ROOT / "paper" / "unitigome_stats.csv")

    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(COLUMNS)
        count = 0
        for row in table_rows(names, mu, means):
            writer.writerow(row)
            count += 1
    print(f"wrote {count} rows to {args.out}")


if __name__ == "__main__":
    main()
