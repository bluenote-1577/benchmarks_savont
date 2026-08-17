#!/usr/bin/env python3
"""Regenerate manuscript-revision figures from compact final result tables.

This plotting-only entry point intentionally does not require raw reads,
alignment files, databases, or completed Savont run directories.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
DEPTHS = [300, 1318, 5792, 25450, 111828]
ABUNDANCE_LABELS = ["<0.02%", "0.02-<0.05%", "0.05-<0.1%", "0.1-<1%", ">=1%"]
DEPTH_LABELS = ["10–<30", "30–<100", "100–<300", "300–<1,000", "≥1,000"]


def set_style(font_size: float = 7, *, compact_axes: bool = False) -> None:
    # Each source figure was originally rendered by a separate process. Reset
    # defaults so generating all panels together does not leak styling between
    # figure groups.
    mpl.rcdefaults()
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": font_size,
            "axes.labelsize": font_size,
            "axes.titlesize": font_size,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    if compact_axes:
        mpl.rcParams.update(
            {
                "axes.linewidth": 0.6,
                "xtick.major.width": 0.6,
                "ytick.major.width": 0.6,
            }
        )


def save_figure(fig: mpl.figure.Figure, prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(
    input_path: Path,
    output_prefix: Path,
    parameter: str,
    parameter_values: list[int] | list[float],
    xlabel: str,
    parameter_labels: list[str] | None = None,
    ylabel: str = "ONT reads",
) -> None:
    results = pd.read_csv(input_path, sep="\t")
    set_style()
    fig, axes = plt.subplots(
        1, 2, figsize=(16.5 / 2.54, 6.7 / 2.54), constrained_layout=True
    )
    for ax, metric, title in zip(
        axes, ["precision", "sensitivity"], ["Precision", "Sensitivity"]
    ):
        matrix = (
            results.pivot(index="read_depth", columns=parameter, values=metric)
            .reindex(index=DEPTHS, columns=parameter_values)
        )
        image = ax.imshow(matrix, vmin=0, vmax=1, cmap="Blues", aspect="auto")
        for row in range(matrix.shape[0]):
            for column in range(matrix.shape[1]):
                value = matrix.iloc[row, column]
                if np.isfinite(value):
                    ax.text(
                        column,
                        row,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        fontsize=5.2,
                        color="white" if value >= 0.62 else "black",
                    )
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_xticks(
            range(len(parameter_values)), parameter_labels or parameter_values
        )
        ax.set_yticks(range(len(DEPTHS)), [f"{depth:,}" for depth in DEPTHS])
        ax.set_ylabel(ylabel)
        for spine in ax.spines.values():
            spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=axes, fraction=0.025, pad=0.02)
    colorbar.set_ticks(np.arange(0, 1.01, 0.2))
    save_figure(fig, output_prefix)


def plot_abundance_fidelity() -> None:
    data = pd.read_csv(
        REPO / "real_notebook/asv_fidelity_by_abundance_subsamples.tsv", sep="\t"
    )
    data["abundance_bin"] = pd.Categorical(
        data["abundance_bin"], categories=ABUNDANCE_LABELS, ordered=True
    )
    means = (
        data.groupby("abundance_bin", observed=False)["pooled_perfect_pct"]
        .mean()
        .reindex(ABUNDANCE_LABELS)
        / 100
    )
    set_style(compact_axes=True)
    fig, ax = plt.subplots(figsize=(8.4 / 2.54, 6.1 / 2.54))
    positions = np.arange(len(ABUNDANCE_LABELS))
    ax.bar(
        positions,
        means,
        width=0.72,
        color="#4C78A8",
        alpha=0.72,
        edgecolor="black",
        linewidth=0.5,
        zorder=1,
    )
    rng = np.random.default_rng(727271)
    for position, abundance_bin in zip(positions, ABUNDANCE_LABELS):
        values = (
            data.loc[data["abundance_bin"] == abundance_bin, "pooled_perfect_pct"]
            .to_numpy()
            / 100
        )
        ax.scatter(
            position + rng.uniform(-0.20, 0.20, size=len(values)),
            values,
            s=12,
            facecolor="#4C78A8",
            edgecolor="black",
            linewidth=0.45,
            alpha=0.95,
            zorder=3,
        )
    ax.set_xticks(positions, ABUNDANCE_LABELS, rotation=35, ha="right")
    ax.set_xlabel("ASV relative abundance within subsample")
    ax.set_ylabel("Fraction of ASVs\nwith a perfect HiFi match")
    ax.set_ylim(0, 1.05)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.45, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout(pad=0.4)
    save_figure(fig, REPO / "real_notebook/asv_fidelity_by_abundance")


def plot_zotu_recall() -> None:
    folder = HERE / "zotu_recall_by_depth"
    datasets = pd.read_csv(folder / "zotu_recall_by_dataset.tsv", sep="\t")
    summary = pd.read_csv(folder / "zotu_recall_by_depth_summary.tsv", sep="\t")
    means = summary.set_index("depth_bin")["mean_dataset_recall"].reindex(DEPTH_LABELS)
    set_style(6, compact_axes=True)
    fig, ax = plt.subplots(figsize=(8.4 / 2.54, 6.3 / 2.54))
    positions = np.arange(len(DEPTH_LABELS))
    ax.bar(
        positions,
        means,
        width=0.72,
        color="#4C78A8",
        alpha=0.72,
        edgecolor="black",
        linewidth=0.5,
        zorder=1,
    )
    rng = np.random.default_rng(727271)
    for position, depth_bin in zip(positions, DEPTH_LABELS):
        values = datasets.loc[
            datasets["depth_bin"].astype(str) == depth_bin, "recall"
        ].to_numpy()
        ax.scatter(
            position + rng.uniform(-0.20, 0.20, size=len(values)),
            values,
            s=12,
            facecolor="#4C78A8",
            edgecolor="black",
            linewidth=0.45,
            alpha=0.95,
            zorder=3,
        )
        mean = means.loc[depth_bin]
        if math.isfinite(mean):
            ax.annotate(
                f"{mean:.2f}",
                (position, mean),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=6,
                bbox={
                    "boxstyle": "round,pad=0.12",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.85,
                },
                zorder=4,
            )
    ax.set_xticks(positions, DEPTH_LABELS, rotation=35, ha="right")
    ax.set_xlabel("Observed ASV depth in ONT sample")
    ax.set_ylabel("Fraction of HiFi UNOISE ASVs\nrecovered with 0 mismatches")
    ax.set_ylim(0, 1.08)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.45, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout(pad=0.4)
    save_figure(fig, folder / "zotu_recall_by_observed_depth")
    # Retain the manuscript's earlier compatibility filename.
    save_figure_copy(
        folder / "zotu_recall_by_observed_depth",
        folder / "zotu_recall_by_expected_depth",
    )


def save_figure_copy(source_prefix: Path, destination_prefix: Path) -> None:
    for suffix in (".svg", ".pdf"):
        destination_prefix.with_suffix(suffix).write_bytes(
            source_prefix.with_suffix(suffix).read_bytes()
        )


def plot_expected_genera() -> None:
    folder = HERE / "its_18s_benchmark"
    data = pd.read_csv(folder / "expected_genera_recovery.tsv", sep="\t")
    methods = ["savont", "unoise", "unoise_min3", "banana"]
    labels = {
        "savont": "Savont",
        "unoise": "UNOISE3",
        "unoise_min3": "UNOISE3 (min. 3 reads)",
        "banana": "BaNaNA\n(OTU only)",
    }
    markers = {"savont": "o", "unoise": "s", "unoise_min3": "^", "banana": "X"}
    linestyles = {"savont": "-", "unoise": "--", "unoise_min3": "-.", "banana": ":"}
    colors = {
        "savont": "#4878d0",
        "unoise": "#ee854a",
        "unoise_min3": "#6acc64",
        "banana": "#956cb4",
    }
    mpl.rcdefaults()
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
    panels = [
        ("18s", "Full-length 18S rRNA protist community\n(ONT R10.4 sup)", 6),
        ("its", "ITS fungal community\n(ONT R10.4 hac)", 10),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(15 / 2.54, 6 / 2.54), sharex=True)
    for ax, (marker, title, expected) in zip(axes, panels):
        subset = data[data["marker"] == marker]
        for method in methods:
            values = subset[subset["method"] == method].sort_values("reads")
            if values.empty:
                continue
            ax.plot(
                values["reads"],
                values["genera_recovered"],
                color=colors[method],
                linestyle=linestyles[method],
                marker=markers[method],
                markersize=4,
                linewidth=1,
                label=labels[method],
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
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.08),
        ncol=5,
        frameon=False,
        fontsize=7,
    )
    fig.tight_layout()
    save_figure(fig, folder / "figures/expected_genera_recovered")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only",
        nargs="+",
        choices=["abundance", "min-cluster", "primary-threshold", "accuracy", "zotu", "its-18s"],
        help="Regenerate only the named figure groups (default: all).",
    )
    args = parser.parse_args()
    selected = set(args.only or [
        "abundance", "min-cluster", "primary-threshold", "accuracy", "zotu", "its-18s"
    ])
    if "abundance" in selected:
        plot_abundance_fidelity()
    if "min-cluster" in selected:
        folder = HERE / "min_cluster_sensitivity"
        plot_heatmap(
            folder / "min_cluster_sensitivity.tsv",
            folder / "min_cluster_sensitivity_heatmaps",
            "min_cluster_size",
            list(range(4, 21, 2)),
            "Minimum cluster size",
        )
    if "primary-threshold" in selected:
        folder = HERE / "primary_clustering_threshold_sensitivity"
        thresholds = [0.91, 0.93, 0.95, 0.97, 0.99]
        plot_heatmap(
            folder / "primary_clustering_threshold_sensitivity.tsv",
            folder / "primary_clustering_threshold_sensitivity_heatmaps",
            "primary_clustering_threshold",
            thresholds,
            "Primary clustering threshold",
            [f"{value:.2f}" for value in thresholds],
        )
    if "accuracy" in selected:
        folder = HERE / "read_accuracy_sensitivity"
        plot_heatmap(
            folder / "read_accuracy_sensitivity.tsv",
            folder / "read_accuracy_sensitivity_heatmaps",
            "mean_read_accuracy_pct",
            [96, 97, 98, 99],
            "Mean simulated read accuracy (%)",
            ylabel="Simulated ONT reads",
        )
    if "zotu" in selected:
        plot_zotu_recall()
    if "its-18s" in selected:
        plot_expected_genera()


if __name__ == "__main__":
    main()
