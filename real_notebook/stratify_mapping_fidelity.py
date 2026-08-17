#!/usr/bin/env python3
"""Stratify ONT Savont-to-HiFi mapping fidelity by ASV abundance.

This uses the same 16 ONT subsamples as the manuscript's pseudo-reference
validation (four environments and four fractions below 0.25). Each ASV call in
each subsample is one observation. ASV abundance is the Savont-estimated depth
encoded in the FASTA header, normalized by the sum of ASV depths in that
subsample.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parents[1]
SAMPLES = ["AD1_ont", "WWTP1_ont", "Zfecal1_ont", "Soil1_ont"]
FRACTIONS = [0.111, 0.037, 0.012, 0.004]

# Quasi-logarithmic relative-abundance bins chosen to resolve the rarest ASVs.
ABUNDANCE_EDGES = [0, 0.0002, 0.0005, 0.001, 0.01, np.inf]
ABUNDANCE_LABELS = [
    "<0.02%",
    "0.02-<0.05%",
    "0.05-<0.1%",
    "0.1-<1%",
    ">=1%",
]

DEPTH_RE = re.compile(r"_depth_(\d+)(?:_|$)")


def set_figure_style() -> None:
    """Match the compact typography and restrained palette of the preprint."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_asv_depths(path: Path) -> dict[str, int]:
    depths: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            asv = line[1:].strip().split()[0]
            match = DEPTH_RE.search(asv)
            if match is None:
                raise ValueError(f"Missing _depth_N in FASTA header: {asv}")
            depths[asv] = int(match.group(1))
    return depths


def read_best_nm(path: Path) -> dict[str, int]:
    """Return the minimum NM tag for each PAF query."""
    best_nm: dict[str, int] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {line_number} in {path}")
            tags = {
                field.split(":", 2)[0]: field.split(":", 2)[2]
                for field in fields[12:]
                if field.count(":") >= 2
            }
            if "NM" not in tags:
                raise ValueError(f"Missing NM tag on PAF line {line_number} in {path}")
            query = fields[0]
            nm = int(tags["NM"])
            best_nm[query] = min(nm, best_nm.get(query, nm))
    return best_nm


def wilson_interval(successes: int, total: int, z: float = 1.95996398454) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    half_width = z * math.sqrt(
        proportion * (1 - proportion) / total + z * z / (4 * total * total)
    ) / denominator
    return center - half_width, center + half_width


def build_observations() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for sample in SAMPLES:
        for fraction in FRACTIONS:
            fasta = BASE / "results_real" / "savont" / f"{sample}_frac{fraction}" / "final_asvs.fasta"
            paf = BASE / "results_real" / "map_ont_pb" / "savont" / f"{sample}_frac{fraction}.map.tsv"
            depths = read_asv_depths(fasta)
            best_nm = read_best_nm(paf)
            total_depth = sum(depths.values())

            unexpected_queries = set(best_nm) - set(depths)
            if unexpected_queries:
                raise ValueError(f"PAF queries absent from {fasta}: {sorted(unexpected_queries)[:3]}")

            for asv, depth in depths.items():
                nm = best_nm.get(asv)
                if nm is None:
                    category = "unmapped"
                elif nm == 0:
                    category = "NM_0"
                elif nm == 1:
                    category = "NM_1"
                else:
                    category = "NM_ge_2"
                rows.append(
                    {
                        "sample": sample,
                        "fraction": fraction,
                        "asv": asv,
                        "asv_depth": depth,
                        "total_asv_depth": total_depth,
                        "relative_abundance": depth / total_depth,
                        "relative_abundance_pct": 100 * depth / total_depth,
                        "NM": nm,
                        "mapping_category": category,
                    }
                )

    observations = pd.DataFrame(rows)
    observations["abundance_bin"] = pd.cut(
        observations["relative_abundance"],
        bins=ABUNDANCE_EDGES,
        labels=ABUNDANCE_LABELS,
        right=False,
        include_lowest=True,
    )
    return observations


def summarize_group(group: pd.DataFrame) -> dict[str, float | int]:
    total = len(group)
    perfect = int((group["mapping_category"] == "NM_0").sum())
    nm_1 = int((group["mapping_category"] == "NM_1").sum())
    nm_ge_2 = int((group["mapping_category"] == "NM_ge_2").sum())
    unmapped = int((group["mapping_category"] == "unmapped").sum())
    ci_low, ci_high = wilson_interval(perfect, total)

    per_subsample = (
        group.assign(perfect=group["mapping_category"].eq("NM_0"))
        .groupby(["sample", "fraction"], observed=True)["perfect"]
        .mean()
    )
    return {
        "asv_observations": total,
        "perfect_n": perfect,
        "NM_1_n": nm_1,
        "NM_ge_2_n": nm_ge_2,
        "unmapped_n": unmapped,
        "pooled_perfect_pct": 100 * perfect / total,
        "pooled_wilson_95_ci_low_pct": 100 * ci_low,
        "pooled_wilson_95_ci_high_pct": 100 * ci_high,
        "subsamples_with_asvs": len(per_subsample),
        "mean_subsample_perfect_pct": 100 * per_subsample.mean(),
        "median_subsample_perfect_pct": 100 * per_subsample.median(),
        "min_subsample_perfect_pct": 100 * per_subsample.min(),
        "max_subsample_perfect_pct": 100 * per_subsample.max(),
    }


def build_summary(observations: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for abundance_bin, group in observations.groupby("abundance_bin", observed=False):
        if group.empty:
            continue
        rows.append({"abundance_bin": str(abundance_bin), **summarize_group(group)})
    rows.append({"abundance_bin": "Overall", **summarize_group(observations)})
    return pd.DataFrame(rows)


def build_subsample_summary(observations: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, fraction, abundance_bin), group in observations.groupby(
        ["sample", "fraction", "abundance_bin"], observed=True
    ):
        row = {
            "sample": sample,
            "fraction": fraction,
            "abundance_bin": str(abundance_bin),
        }
        row.update(summarize_group(group))
        rows.append(row)
    return pd.DataFrame(rows)


def plot_subsample_fidelity(subsample_summary: pd.DataFrame, output_prefix: Path) -> None:
    """Bar means with one dot for every sample-by-depth observation."""
    set_figure_style()
    plot_data = subsample_summary.copy()
    plot_data["abundance_bin"] = pd.Categorical(
        plot_data["abundance_bin"], categories=ABUNDANCE_LABELS, ordered=True
    )

    means = (
        plot_data.groupby("abundance_bin", observed=False)["pooled_perfect_pct"]
        .mean()
        .reindex(ABUNDANCE_LABELS)
        / 100
    )

    fig, ax = plt.subplots(figsize=(8.4 / 2.54, 6.1 / 2.54))
    x_positions = np.arange(len(ABUNDANCE_LABELS))
    ax.bar(
        x_positions,
        means,
        width=0.72,
        color="#4C78A8",
        alpha=0.72,
        edgecolor="black",
        linewidth=0.5,
        zorder=1,
    )

    rng = np.random.default_rng(727271)
    for x_position, abundance_bin in zip(x_positions, ABUNDANCE_LABELS):
        values = (
            plot_data.loc[
                plot_data["abundance_bin"] == abundance_bin, "pooled_perfect_pct"
            ].to_numpy()
            / 100
        )
        jitter = rng.uniform(-0.20, 0.20, size=len(values))
        ax.scatter(
            x_position + jitter,
            values,
            s=12,
            facecolor="#4C78A8",
            edgecolor="black",
            linewidth=0.45,
            alpha=0.95,
            zorder=3,
        )

    ax.set_xticks(x_positions, ABUNDANCE_LABELS, rotation=35, ha="right")
    ax.set_xlabel("ASV relative abundance within subsample")
    ax.set_ylabel("Fraction of ASVs\nwith a perfect HiFi match")
    ax.set_ylim(0, 1.05)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.45, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout(pad=0.4)
    fig.savefig(output_prefix.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    observations = build_observations()
    summary = build_summary(observations)
    subsample_summary = build_subsample_summary(observations)

    output_prefix = Path(__file__).resolve().parent / "asv_fidelity_by_abundance"
    observations.to_csv(f"{output_prefix}_observations.tsv", sep="\t", index=False)
    summary.to_csv(f"{output_prefix}_summary.tsv", sep="\t", index=False, float_format="%.4f")
    subsample_summary.to_csv(
        f"{output_prefix}_subsamples.tsv", sep="\t", index=False, float_format="%.4f"
    )
    plot_subsample_fidelity(subsample_summary, output_prefix)

    print(summary.to_string(index=False, float_format=lambda value: f"{value:.2f}"))
    print(
        f"\nWrote {len(observations):,} ASV observations from "
        f"{len(SAMPLES) * len(FRACTIONS)} subsamples, plus SVG and PDF figures."
    )


if __name__ == "__main__":
    main()
