import argparse
import sys
from matplotlib.ticker import FuncFormatter, LogFormatterSciNotation, LogLocator, NullFormatter
import numpy as np

from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
ROOT = TEST_DIR.parent
sys.path.insert(0, str(ROOT))

GENOME_LENGTH = 13.4e6  # 13.4 M is genome length

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
    rho = B/(gamma*L*GENOME_LENGTH)
    m1_z_list = [get_z_val(1 / kmer_means[k], k, rho, L) for k in kmer_list]
    return max(m1_z_list)

def get_m2_z_rho(B, kmer_means, kmer_list, L = 125):
    rho = B/(L*GENOME_LENGTH)
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


def get_crossings(x_list, list1, list2):
    """(x, z) where the two curves cross, interpolated linearly in log-log space."""
    lx, d = np.log(x_list), np.log(list1) - np.log(list2)
    crossings = []
    for i in get_sign_change_positions(list1, list2):  # sign changes between i-1 and i
        t = d[i - 1] / (d[i - 1] - d[i])
        x = np.exp(lx[i - 1] + t * (lx[i] - lx[i - 1]))
        z = np.exp(np.log(list1[i - 1]) + t * (np.log(list1[i]) - np.log(list1[i - 1])))
        crossings.append((x, z))
    return crossings


def add_coverage_axis(ax, location, cov_per_x, label):
    """Extra x axis at location ('top', 'bottom' or an offset) showing coverage = x * cov_per_x."""
    sec = ax.secondary_xaxis(location, functions=(lambda x: x * cov_per_x, lambda c: c / cov_per_x))
    sec.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    sec.set_xlabel(label, fontsize=15, labelpad=10)  # gap between the label and its tick labels
    return sec


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
    fig,axs = plt.subplots(1, 3,figsize=(12,6.6), sharey=True, sharex=True,
                             squeeze=False)
    axs = axs.flatten()
    m2_cov_list = budget_list / GENOME_LENGTH  # x axis: M2 coverage = B / G
    gamma_list = [1, 5, 10]
    for i, (ax, gamma) in enumerate(zip(axs, gamma_list)):
        m1_z_list = [get_m1_z_rho(B, kmer_means, m1_kmer_list, gamma, L1) for B in budget_list]
        m2_z_list = [get_m2_z_rho(B, kmer_means, m2_kmer_list, L2) for B in budget_list]
        crossings = get_crossings(m2_cov_list, m1_z_list, m2_z_list)
        if not crossings:
            print(f"gamma = {gamma:g}: no crossing")
        for x, z in crossings:
            # x is M2 coverage (B / G); M1 coverage is B / (gamma * G) = x / gamma
            print(f"gamma = {gamma:g}: M2 coverage = {x:.2f}, M1 coverage = {x / gamma:.2f}, "
                  f"Z = {z:,.0f}")
            # through the crossing, from the bottom axis (M2 coverage) to the top axis (M1 coverage)
            ax.axvline(x, color="#CC79A7", linestyle="--", lw=1.5)
            ax.plot(x, z, "o", ms=8, color="#CC79A7", mec="white", mew=1.5, zorder=5)
        ax.plot(m2_cov_list, m1_z_list, ls='-', lw=2.5, color="#0072B2", label=rf"$M_1$ ($L$ = {L1})")
        ax.plot(m2_cov_list, m2_z_list, ls='-', lw=2.5, color="#E69F00", label=rf"$M_2$ ($L$ = {L2})")
        ax.text(0.97, 0.05, "ABCDEF"[i], transform=ax.transAxes, fontsize=22,
                    fontweight="bold", ha="right", va="bottom", color="#222")
        ax.set_xlabel(r"$M_2$ coverage", fontsize=15)
        # ax.spines["bottom"].set_color("#E69F00")
        # ax.spines["bottom"].set_linewidth(2)
        ax.set_xscale('log')
        add_coverage_axis(ax, "top", 1 / gamma, r"$M_1$ coverage")
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(LogLocator(subs=(1, 2, 5)))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.grid(True, which="both", alpha=0.25)
        ax.set_title(rf"$\gamma = {gamma:g}$", pad=24)  # extra room above the M1 coverage label
    axs[0].set_ylabel(r"$Z_\rho$-value")
    
    fig.subplots_adjust(left=0.11, right=0.98, top=0.78, bottom=0.19, wspace=0.08)
    handles, labels = axs[0].get_legend_handles_labels()
    x_mid = (axs[0].get_position().x0 + axs[-1].get_position().x1) / 2
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(x_mid, 0), ncol=2, frameon=False)
    fig.savefig(out_fn)
    print(f"Plot saved to {out_fn}")



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", default=ROOT / "paper" / "p_brunnea_reftig_stats.txt")
    ap.add_argument("--out", default=ROOT / "paper" /"figure_5.pdf")
    args = ap.parse_args()
    fn = args.stats
    out_fn = args.out
    kmer_means = get_kmer_means(fn)
    
    budget_list = np.logspace(0, 3, 61)  
    budget_list *= GENOME_LENGTH
    
    plot_z_vs_B(budget_list, kmer_means, out_fn)

main()