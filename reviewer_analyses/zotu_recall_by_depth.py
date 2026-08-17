#!/usr/bin/env python3
"""Estimate Savont recall across directly observed HiFi-zOTU depth bins.

For every ONT subsample, the trimmed reads are dereplicated only to preserve
read multiplicities efficiently and are then mapped directly to that sample's
HiFi + UNOISE3 pseudo-reference with USEARCH ``-otutab``. No ONT-derived
UNOISE3 zOTUs are generated or used. A pseudo-reference zOTU is recovered when
at least one Savont ASV has an NM=0 alignment to it. Consequently, the depth
used for binning is the realized mapped-read count in the same ONT subsample,
not an expectation obtained by scaling a full-depth count.
"""

from __future__ import annotations

import math
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "zotu_recall_by_depth"
USEARCH = Path("/homes9/jshaw/bin/usearch")
SAMPLES = {
    "AD1_ont": "AD",
    "Soil1_ont": "Soil",
    "WWTP1_ont": "WWTP",
    "Zfecal1_ont": "Zfecal",
}
FRACTIONS = [1.0, 0.33, 0.111, 0.037, 0.012, 0.004]
DEPTH_EDGES = [10, 30, 100, 300, 1000, np.inf]
DEPTH_LABELS = ["10–<30", "30–<100", "100–<300", "300–<1,000", "≥1,000"]
PARALLEL_DATASETS = 4


def read_otutab(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    table = table.rename(columns={table.columns[0]: "zotu"})
    count_columns = [column for column in table.columns if column != "zotu"]
    table["observed_depth"] = table[count_columns].sum(axis=1)
    return table[["zotu", "observed_depth"]]


def count_fastq(path: Path) -> int:
    with path.open() as handle:
        lines = sum(1 for _ in handle)
    if lines % 4:
        raise ValueError(f"Incomplete FASTQ record in {path}")
    return lines // 4


def direct_hifi_otutab(sample: str, fraction: float) -> tuple[pd.DataFrame, dict[str, object]]:
    """Map one ONT subsample directly to its HiFi pseudo-reference."""
    dataset = f"{sample}_frac{fraction:g}"
    reads = BASE / "results_real" / "subsampled" / f"{dataset}.fastq"
    reference_name = sample.replace("1_ont", "")
    reference = BASE / "results_real" / "pb_cat_refs" / f"{reference_name}_pb_zotus.fa"
    direct_dir = OUT / "direct_hifi_otutabs"
    log_dir = OUT / "logs"
    direct_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    uniques = direct_dir / f"{dataset}.uniques.fa"
    otutab_path = direct_dir / f"{dataset}.tsv"
    derep_log = log_dir / f"{dataset}.dereplicate.log"
    otutab_log = log_dir / f"{dataset}.otutab.log"

    # The OTU table is the durable result. If it exists, resume without
    # regenerating the large, losslessly dereplicated intermediate FASTA.
    if not otutab_path.exists():
        with derep_log.open("w") as log_handle:
            subprocess.run(
                [
                    str(USEARCH),
                    "-fastx_uniques",
                    str(reads),
                    "-sizeout",
                    "-relabel",
                    "Read_",
                    "-minuniquesize",
                    "1",
                    "-fastaout",
                    str(uniques),
                ],
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                check=True,
            )
        with otutab_log.open("w") as log_handle:
            subprocess.run(
                [
                    str(USEARCH),
                    "-otutab",
                    str(uniques),
                    "-zotus",
                    str(reference),
                    "-otutabout",
                    str(otutab_path),
                    "-sample_delim",
                    "_",
                    "-threads",
                    "20",
                ],
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                check=True,
            )

    depths = read_otutab(otutab_path)
    input_reads = count_fastq(reads)
    mapped_reads = int(depths["observed_depth"].sum())
    depths.insert(0, "dataset", dataset)
    depths.insert(0, "fraction", fraction)
    depths.insert(0, "sample", sample)
    audit = {
        "sample": sample,
        "fraction": fraction,
        "dataset": dataset,
        "input_reads": input_reads,
        "mapped_reads": mapped_reads,
        "mapped_fraction": mapped_reads / input_reads,
        "hifi_zotus_with_mapped_reads": len(depths),
    }
    return depths, audit


def read_recovered_targets(path: Path) -> set[str]:
    recovered: set[str] = set()
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            nm = next(
                (int(tag.rsplit(":", 1)[1]) for tag in fields[12:] if tag.startswith("NM:i:")),
                None,
            )
            if nm == 0:
                recovered.add(fields[5])
    return recovered


def build_observations(depths: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for sample in SAMPLES:
        for fraction in FRACTIONS:
            dataset = f"{sample}_frac{fraction:g}"
            dataset_depths = depths.loc[depths["dataset"] == dataset]
            paf = (
                BASE
                / "results_real"
                / "map_ont_pb"
                / "savont"
                / f"{sample}_frac{fraction:g}.map.tsv"
            )
            recovered = read_recovered_targets(paf)
            for record in dataset_depths.itertuples(index=False):
                rows.append(
                    {
                        "sample": sample,
                        "fraction": fraction,
                        "dataset": dataset,
                        "zotu": record.zotu,
                        "observed_depth": record.observed_depth,
                        "recovered_NM0": record.zotu in recovered,
                    }
                )

    observations = pd.DataFrame(rows)
    observations["depth_bin"] = pd.cut(
        observations["observed_depth"],
        bins=DEPTH_EDGES,
        labels=DEPTH_LABELS,
        right=False,
    )
    return observations


def summarize_datasets(observations: pd.DataFrame) -> pd.DataFrame:
    in_bins = observations.dropna(subset=["depth_bin"])
    summary = (
        in_bins.groupby(["sample", "fraction", "dataset", "depth_bin"], observed=True)
        .agg(zotus=("zotu", "size"), recovered_NM0=("recovered_NM0", "sum"))
        .reset_index()
    )
    summary["recall"] = summary["recovered_NM0"] / summary["zotus"]
    return summary


def summarize_bins(dataset_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for depth_bin, group in dataset_summary.groupby("depth_bin", observed=False):
        if group.empty:
            continue
        rows.append(
            {
                "depth_bin": str(depth_bin),
                "dataset_points": len(group),
                "mean_dataset_recall": group["recall"].mean(),
                "sd_dataset_recall": group["recall"].std(ddof=1),
                "median_dataset_recall": group["recall"].median(),
                "min_dataset_recall": group["recall"].min(),
                "max_dataset_recall": group["recall"].max(),
                "zotu_observations": int(group["zotus"].sum()),
                "recovered_NM0": int(group["recovered_NM0"].sum()),
                "pooled_recall": group["recovered_NM0"].sum() / group["zotus"].sum(),
            }
        )
    return pd.DataFrame(rows)


def set_figure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 6,
            "axes.labelsize": 6,
            "axes.titlesize": 6,
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


def plot_recall(dataset_summary: pd.DataFrame, bin_summary: pd.DataFrame) -> None:
    set_figure_style()
    fig, ax = plt.subplots(figsize=(8.4 / 2.54, 6.3 / 2.54))
    means = bin_summary.set_index("depth_bin")["mean_dataset_recall"].reindex(DEPTH_LABELS)
    x_positions = np.arange(len(DEPTH_LABELS))
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
    for x_position, depth_bin in zip(x_positions, DEPTH_LABELS):
        values = dataset_summary.loc[
            dataset_summary["depth_bin"].astype(str) == depth_bin, "recall"
        ].to_numpy()
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
        mean = means.loc[depth_bin]
        if math.isfinite(mean):
            ax.annotate(
                f"{mean:.2f}",
                (x_position, mean),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=6,
                bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.85},
                zorder=4,
            )

    ax.set_xticks(x_positions, DEPTH_LABELS, rotation=35, ha="right")
    ax.set_xlabel("Observed ASV depth in ONT sample")
    ax.set_ylabel("Fraction of HiFi UNOISE ASVs\nrecovered with 0 mismatches")
    ax.set_ylim(0, 1.08)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.45, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout(pad=0.4)
    fig.savefig(OUT / "zotu_recall_by_observed_depth.svg", bbox_inches="tight")
    fig.savefig(OUT / "zotu_recall_by_observed_depth.pdf", bbox_inches="tight")
    # Compatibility outputs for the manuscript source while it is being revised.
    fig.savefig(OUT / "zotu_recall_by_expected_depth.svg", bbox_inches="tight")
    fig.savefig(OUT / "zotu_recall_by_expected_depth.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    depth_tables = []
    audits = []
    datasets = [(sample, fraction) for sample in SAMPLES for fraction in FRACTIONS]
    with ThreadPoolExecutor(max_workers=PARALLEL_DATASETS) as executor:
        futures = {
            executor.submit(direct_hifi_otutab, sample, fraction): (sample, fraction)
            for sample, fraction in datasets
        }
        for future in as_completed(futures):
            sample, fraction = futures[future]
            depths, audit = future.result()
            depth_tables.append(depths)
            audits.append(audit)
            print(
                f"Completed {sample}_frac{fraction:g}: "
                f"{audit['mapped_reads']:,}/{audit['input_reads']:,} reads assigned "
                f"({audit['mapped_fraction']:.1%})",
                flush=True,
            )
    depth_table = (
        pd.concat(depth_tables, ignore_index=True)
        .sort_values(["sample", "fraction", "zotu"], ascending=[True, False, True])
        .reset_index(drop=True)
    )
    audits = sorted(audits, key=lambda row: (str(row["sample"]), -float(row["fraction"])))
    observations = build_observations(depth_table)
    dataset_summary = summarize_datasets(observations)
    bin_summary = summarize_bins(dataset_summary)
    depth_table.to_csv(OUT / "hifi_zotu_observed_depths.tsv", sep="\t", index=False)
    pd.DataFrame(audits).to_csv(OUT / "direct_mapping_audit.tsv", sep="\t", index=False)
    observations.to_csv(OUT / "zotu_recall_observations.tsv", sep="\t", index=False)
    dataset_summary.to_csv(
        OUT / "zotu_recall_by_dataset.tsv", sep="\t", index=False, float_format="%.6f"
    )
    bin_summary.to_csv(
        OUT / "zotu_recall_by_depth_summary.tsv", sep="\t", index=False, float_format="%.6f"
    )
    plot_recall(dataset_summary, bin_summary)
    print(bin_summary.to_string(index=False))


if __name__ == "__main__":
    main()
