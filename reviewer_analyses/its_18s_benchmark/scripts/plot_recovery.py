#!/usr/bin/env python3
"""Plot expected-genus recovery (supplementary figure), styled to match the
manuscript's Figure 2.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# BaNaNA results remain in the source tables for supplementary reporting but
# are not included in this figure.
METHODS = ["savont", "unoise", "unoise_min3"]
LABELS = {
    "savont": "Savont",
    "unoise": "UNOISE3",
    "unoise_min3": "UNOISE3 (min. 3 reads)",
}
MARKERS = {"savont": "o", "unoise": "s", "unoise_min3": "^"}
LINESTYLES = {"savont": "-", "unoise": "--", "unoise_min3": "-."}


def set_style() -> dict[str, tuple[float, float, float]]:
    palette = sns.color_palette("muted")
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 7,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    # palette[3] is the red used for DADA2 in Figure 2; skipped here since
    # DADA2 is not run on ONT ITS/18S data.
    return {
        "savont": palette[0],
        "unoise": palette[1],
        "unoise_min3": palette[2],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    data = pd.read_csv(args.input, sep="\t")
    colors = set_style()
    panels = [
        ("18s", "Full-length 18S rRNA protist community\n(ONT R10.4 sup)", 6),
        ("its", "ITS fungal community\n(ONT R10.4 hac)", 10),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(15 / 2.54, 6 / 2.54), sharex=True)
    for ax, (marker, title, expected) in zip(axes, panels):
        subset = data[data["marker"] == marker]
        for method in METHODS:
            values = subset[subset["method"] == method].sort_values("reads")
            if values.empty:
                continue
            ax.plot(
                values["reads"],
                values["genera_recovered"],
                color=colors[method],
                linestyle=LINESTYLES[method],
                marker=MARKERS[method],
                markersize=4,
                linewidth=1,
                label=LABELS[method],
            )
        ax.set_xscale("log")
        ax.set_title(title)
        ax.set_xlabel("# reads")
        ax.set_ylabel("Expected genera recovered")
        ax.set_ylim(-0.25, expected + 0.45)
        ax.set_yticks(range(0, expected + 1))
        ax.axhline(
            expected,
            color="0.35",
            linestyle="--",
            linewidth=1,
            label="Expected # of genera in sample",
        )
        ax.grid(True, alpha=0.3, linestyle="--")
        ax.set_axisbelow(True)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.08),
        ncol=4,
        frameon=False,
        fontsize=7,
    )
    fig.tight_layout()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(args.output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
