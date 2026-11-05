import argparse
import csv
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

TEST_DIR = Path(__file__).resolve().parent
ROOT = TEST_DIR.parent
sys.path.insert(0, str(ROOT))

DATA = ROOT / "paper" / "figure_6_data"
SOURCES = {  
    "bcalm": (DATA / "bcalm_reftig_stats.txt", DATA / "unitigome_stat.csv",
              "#0072B2", "o"), 
    "Velvet": (DATA / "vetvet_reftig_stat.txt", DATA / "velvet_contig_stats.csv",
               "#D55E00", "D"),
}
SPECIES = {"GCF_037157775.1": "P. brunnea", "GCF_060026735.1": "J. thaxteri",
           "GCF_009914755.1": "H. sapiens (Chr 20)"} 
ITALIC = {"GCF_037157775.1", "GCF_060026735.1",  "GCF_009914755.1"} 
OBS_LABEL = {"bcalm": "Alg. 1", "Velvet": "Velvet"}  
CURVE_LABEL = {"bcalm": r"$\mu_k$", "Velvet": r"$\tilde \mu_k$"}
INK = "#333333"



plt.rcParams.update({"font.size": 17, "axes.labelsize": 18, "xtick.labelsize": 15,
                     "ytick.labelsize": 15, "legend.fontsize": 14,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

def expected_z(mu, rho, L, k):
    ell = L - k + 1
    num = (rho + mu) * np.exp(-mu * L - rho * ell)
    den = mu + rho * np.exp(-(rho + mu) * ell)
    return (1.0 - num / den) / mu

def accession_of(name):
    return "_".join(name.replace(" ", "_").split("_")[:2])


def load_means(path, accession):
    means = {}
    with open(path) as f:
        for line in f:
            fields = line.split(maxsplit=4)
            if len(fields) >= 4 and accession_of(fields[0]) == accession:
                means[int(fields[1])] = float(fields[3])
    return means


def predict(means, ks, covs, L):
    ez = np.array([[expected_z(1.0 / means[k], c / L, L, k) for k in ks] for c in covs])
    return ez.max(axis=1), np.array(ks)[ez.argmax(axis=1)]


def observed_maxima(path, accession, ks, col="Mean"):
    best = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            k = int(row["k-mer"]) if "k-mer" in row else None
            if accession_of(row["Genome"]) != accession or (k is not None and k not in ks):
                continue
            cov, mean = float(row["Coverage"]), float(row[col])
            if cov not in best or mean > best[cov][0]:
                best[cov] = (mean, k)
    return dict(sorted(best.items()))


def draw_panel(ax, acc, inputs, L, kmin, kmax, covs):
    rows, handles, curves, limits = [], [], [], []
    for name, (mu_path, obs_path) in inputs.items():
        color, marker = SOURCES[name][2:]
        means = load_means(mu_path, acc)
        ks = [k for k in sorted(means) if kmin <= k <= min(kmax, L - 1)]
        pred, _ = predict(means, ks, covs, L)
        limit = means[ks[-1]]
        limits.append(limit)
        curves += ax.plot(covs, pred, "-", lw=2.5, color=color, label=CURVE_LABEL[name])
        obs = observed_maxima(obs_path, acc, set(ks))
        obs_c = np.array(list(obs))
        handles += ax.plot(obs_c, [m for m, _ in obs.values()], marker, ms=9, mfc="none",
                           mew=2, color=color, zorder=4, label=OBS_LABEL[name])
        p_obs, k_obs = predict(means, ks, obs_c, L)
        rows += [(acc, name, c, m, k, p, kp) for c, (m, k), p, kp in zip(obs_c, obs.values(), p_obs, k_obs)]
    ax.set_xlabel("Coverage")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_xlim(0, 260)
    ax.grid(True, alpha=0.25)
    return handles, curves, rows, 1.15 * max(limits)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bcalm-reftigs", default=SOURCES["bcalm"][0])
    ap.add_argument("--bcalm-unitigs", default=SOURCES["bcalm"][1])
    ap.add_argument("--velvet-reftigs", default=SOURCES["Velvet"][0])
    ap.add_argument("--velvet-contigs", default=SOURCES["Velvet"][1])
    ap.add_argument("--accessions", nargs="+", default=list(SPECIES),
                    help="one panel per genome, left to right")
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--kmin", type=int, default=20)
    ap.add_argument("--kmax", type=int, default=141)
    ap.add_argument("--out", default=ROOT / "paper")
    args = ap.parse_args()

    accs = [accession_of(a) for a in args.accessions]
    L = args.read_length
    inputs = {"bcalm": (args.bcalm_reftigs, args.bcalm_unitigs),
              "Velvet": (args.velvet_reftigs, args.velvet_contigs)}
    covs = np.linspace(2, 250, 1000)
    fig, axes = plt.subplots(1, len(accs), figsize=(16, 6),
                             squeeze=False)
    # rows = []
    for i, (ax, acc) in enumerate(zip(axes[0], accs)):
        handles, curves, _, top = draw_panel(ax, acc, inputs, L, args.kmin, args.kmax, covs)
        # rows += r # if we want to see the data later
        ax.set_ylim(0, top)
        
        ax.set_title(SPECIES.get(acc, acc), style="italic" if acc in ITALIC else "normal")
        ax.text(0.04, 0.95, "ABCDEF"[i], transform=ax.transAxes, fontsize=22,
                fontweight="bold", ha="left", va="top", color=INK)
    axes[0, 0].set_ylabel(r"$Z_\rho$-value")
    fig.legend(handles=handles + curves, loc="lower center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.09, 1, 1))

    out = Path(args.out)
    # out.mkdir(parents=True, exist_ok=True)
    fn = f"figure_6"
    fig.savefig(out / f"{fn}.pdf", dpi=300)
    print(f"Saved figure to {out / f'{fn}.pdf'}")
    

if __name__ == "__main__":
    main()
