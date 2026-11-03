import argparse
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import numpy as np

import sys
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
ROOT = TEST_DIR.parent
sys.path.insert(0, str(ROOT))

def get_kmer_means(fn):
    means = {}
    with open(fn) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if i == 0:
            continue
        fields = line.split(maxsplit=4)
        if len(fields) >= 4:
            means[int(fields[1])] = float(fields[3])
    return means


def get_z_val(mu, k, rho, L):
    ell = L - k + 1
    num = (rho + mu) * np.exp(-mu * L - rho * ell)
    den = mu + rho * np.exp(-(rho + mu) * ell)
    return (1.0 - num / den) / mu

def get_m1_z_rho(B, kmer_means, kmer_list, gamma, L = 300):
    rho = B/(gamma*L)
    m1_z_list = [get_z_val(1 / kmer_means[k], k, rho, L) for k in kmer_list]
    return max(m1_z_list)

def get_m2_z_rho(B, kmer_means, kmer_list, L = 125):
    rho = B/L
    m2_z_list = [get_z_val(1 / kmer_means[k], k, rho, L) for k in kmer_list]
    return max(m2_z_list)

def get_sign_change_positions(list1, list2):
    positions = []
    prev_val = None
    for i, (val1, val2) in enumerate(zip(list1, list2)):
        val = val1 - val2
        if prev_val is not None and val * prev_val < 0:
            positions.append(i)
        prev_val = val
    return positions


def plot_z_vs_B(budget_list, kmer_means, out_fn):
    kmer_list = list(kmer_means.keys())
    L1 = 300
    L2 = 125
    m1_kmer_list = [k for k in kmer_list if k <= L1]
    m2_kmer_list = [k for k in kmer_list if k <= L2]
    
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 17, "axes.labelsize": 18, "xtick.labelsize": 15,
                     "ytick.labelsize": 15, "legend.fontsize": 16,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
    fig,axs = plt.subplots(1, 3,figsize=(12,5), sharey=True, sharex=True,
                             squeeze=False)
    axs = axs.flatten()
    scale = 13.4e6 / 1e9 # Per Giga BP (13.4 M is genome length)
    scaled_budget_list = budget_list * scale
    gamma_list = [1, 2, 3]
    for i, (ax, gamma) in enumerate(zip(axs, gamma_list)):
        m1_z_list = [get_m1_z_rho(B, kmer_means, m1_kmer_list, gamma, L1) for B in budget_list]
        m2_z_list = [get_m2_z_rho(B, kmer_means, m2_kmer_list, L2) for B in budget_list]
        sign_change_positions = get_sign_change_positions(m1_z_list, m2_z_list)
        for pos in sign_change_positions:
            ax.axvline(scaled_budget_list[pos], color="#CC79A7", linestyle='--')
        ax.plot(scaled_budget_list, m1_z_list, ls='-', lw=2.5, color="#0072B2", label=rf"$M_1$ ($L$ = {L1})")
        ax.plot(scaled_budget_list, m2_z_list, ls='-', lw=2.5, color="#E69F00", label=rf"$M_2$ ($L$ = {L2})")
        ax.text(0.97, 0.05, "ABCDEF"[i], transform=ax.transAxes, fontsize=22,
                    fontweight="bold", ha="right", va="bottom", color="#222")
        ax.set_xlabel('Budget')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(LogLocator(subs=(1, 2, 5)))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.grid(True, which="both", alpha=0.25)
        ax.set_title(rf"$\gamma = {gamma:g}$")
    axs[0].set_ylabel(r"$Z_\rho$-value")
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(out_fn)
    print(f"Plot saved to {out_fn}")



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", type=Path, default=TEST_DIR /"p_brunnea_reftig_stats.txt")
    ap.add_argument("--out", type=Path, default=TEST_DIR / "two-machine-budget.pdf")
    args = ap.parse_args()
    fn = args.stats
    out_fn = args.out
    kmer_means = get_kmer_means(fn)
    
    budget_list = np.logspace(0, 3, 61)  
    
    plot_z_vs_B(budget_list, kmer_means, out_fn)

main()