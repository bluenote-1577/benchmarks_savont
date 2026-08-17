#!/usr/bin/env python3
"""Sensitivity analysis for Savont's primary clustering threshold on ONT Zymo data.

The five read depths and exact-match scoring mirror min_cluster_sensitivity.py.
Existing successful runs are reused, so the script can be safely resumed.
"""

from __future__ import annotations

import gzip
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "primary_clustering_threshold_sensitivity"
SAVONT = Path("/homes9/jshaw/.cargo/bin/savont")
MINIMAP2 = Path("/homes9/jshaw/bin/minimap2")
TRUTH = BASE / "morten_data" / "V1V8_bacterial_amplicons.fasta.gz"
PRIMARY_CLUSTERING_THRESHOLDS = [0.91, 0.93, 0.95, 0.97, 0.99]
READ_DEPTH_FILES = {
    300: "v1v8_zmock1_frac300.0000000000001.fastq",
    1318: "v1v8_zmock1_frac1318.1911682282378.fastq",
    5792: "v1v8_zmock1_frac5792.093186649751.fastq",
    25450: "v1v8_zmock1_frac25450.286947322187.fastq",
    111828: "v1v8_zmock1_frac111827.81160944818.fastq",
}
THREADS_PER_RUN = 10
PARALLEL_RUNS = 2


def threshold_label(threshold: float) -> str:
    return f"{threshold:.2f}"


def count_fasta(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return sum(line.startswith(">") for line in handle)


def is_complete(run_dir: Path) -> bool:
    fasta = run_dir / "final_asvs.fasta"
    if not fasta.exists():
        return False
    logs = sorted(run_dir.glob("savont_*.log"))
    return bool(logs) and "COMPLETED SUCCESSFULLY" in logs[-1].read_text(errors="replace")


def run_one(read_depth: int, threshold: float) -> tuple[int, float]:
    reads = BASE / "results-FINAL-preprint" / "subsampled" / READ_DEPTH_FILES[read_depth]
    run_dir = OUT / "results" / f"reads_{read_depth}" / f"threshold_{threshold_label(threshold)}"
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    if is_complete(run_dir):
        return read_depth, threshold
    if run_dir.exists() and any(run_dir.iterdir()):
        raise RuntimeError(f"Incomplete non-empty run directory requires inspection: {run_dir}")

    command = [
        str(SAVONT),
        "asv",
        "-t",
        str(THREADS_PER_RUN),
        str(reads),
        "-o",
        str(run_dir),
        "--primary-clustering-threshold",
        threshold_label(threshold),
    ]
    subprocess.run(command, check=True)
    if not is_complete(run_dir):
        raise RuntimeError(f"Savont exited without a successful completion marker: {run_dir}")
    return read_depth, threshold


def parse_nm(fields: list[str]) -> int:
    for tag in fields[12:]:
        if tag.startswith("NM:i:"):
            return int(tag.rsplit(":", 1)[1])
    raise ValueError("PAF record is missing NM:i")


def score_one(read_depth: int, threshold: float, truth_count: int) -> dict[str, object]:
    label = threshold_label(threshold)
    run_dir = OUT / "results" / f"reads_{read_depth}" / f"threshold_{label}"
    asvs = run_dir / "final_asvs.fasta"
    paf = OUT / "alignments" / f"reads_{read_depth}_threshold_{label}.paf"
    paf.parent.mkdir(parents=True, exist_ok=True)
    command = [str(MINIMAP2), "--cs", "--secondary=no", str(TRUTH), str(asvs)]
    with paf.open("w") as output:
        subprocess.run(command, stdout=output, check=True)

    best_query_nm: dict[str, int] = {}
    exact_targets: set[str] = set()
    with paf.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            nm = parse_nm(fields)
            query, target = fields[0], fields[5]
            best_query_nm[query] = min(nm, best_query_nm.get(query, nm))
            if nm == 0:
                exact_targets.add(target)

    asv_count = count_fasta(asvs)
    exact_asv_count = sum(nm == 0 for nm in best_query_nm.values())
    false_positive_count = asv_count - exact_asv_count
    precision = exact_asv_count / asv_count if asv_count else np.nan
    sensitivity = len(exact_targets) / truth_count
    return {
        "read_depth": read_depth,
        "primary_clustering_threshold": threshold,
        "asv_count": asv_count,
        "exact_asv_count": exact_asv_count,
        "false_positive_count": false_positive_count,
        "truth_asvs_recovered": len(exact_targets),
        "truth_asv_count": truth_count,
        "precision": precision,
        "sensitivity": sensitivity,
    }


def set_figure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_heatmaps(results: pd.DataFrame) -> None:
    set_figure_style()
    fig, axes = plt.subplots(1, 2, figsize=(16.5 / 2.54, 6.7 / 2.54), constrained_layout=True)
    depths = list(READ_DEPTH_FILES)
    for ax, metric, title in zip(
        axes,
        ["precision", "sensitivity"],
        ["Precision", "Sensitivity"],
    ):
        matrix = (
            results.pivot(
                index="read_depth", columns="primary_clustering_threshold", values=metric
            ).reindex(index=depths, columns=PRIMARY_CLUSTERING_THRESHOLDS)
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
        ax.set_xlabel("Primary clustering threshold")
        ax.set_xticks(
            range(len(PRIMARY_CLUSTERING_THRESHOLDS)),
            [threshold_label(value) for value in PRIMARY_CLUSTERING_THRESHOLDS],
        )
        ax.set_yticks(range(len(depths)), [f"{depth:,}" for depth in depths])
        ax.set_ylabel("ONT reads")
        for spine in ax.spines.values():
            spine.set_visible(False)

    colorbar = fig.colorbar(image, ax=axes, fraction=0.025, pad=0.02)
    colorbar.set_ticks(np.arange(0, 1.01, 0.2))
    fig.savefig(OUT / "primary_clustering_threshold_sensitivity_heatmaps.svg", bbox_inches="tight")
    fig.savefig(OUT / "primary_clustering_threshold_sensitivity_heatmaps.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    version = subprocess.run([str(SAVONT), "--version"], check=True, capture_output=True, text=True)
    (OUT / "software_version.txt").write_text(version.stdout)

    jobs = [
        (depth, threshold)
        for depth in READ_DEPTH_FILES
        for threshold in PRIMARY_CLUSTERING_THRESHOLDS
    ]
    with ThreadPoolExecutor(max_workers=PARALLEL_RUNS) as executor:
        futures = {executor.submit(run_one, *job): job for job in jobs}
        for future in as_completed(futures):
            depth, threshold = future.result()
            print(
                f"Completed/read existing: reads={depth:,}, threshold={threshold_label(threshold)}",
                flush=True,
            )

    truth_count = count_fasta(TRUTH)
    rows = [score_one(depth, threshold, truth_count) for depth, threshold in jobs]
    results = pd.DataFrame(rows).sort_values(["read_depth", "primary_clustering_threshold"])
    results.to_csv(
        OUT / "primary_clustering_threshold_sensitivity.tsv",
        sep="\t",
        index=False,
        float_format="%.6f",
    )
    plot_heatmaps(results)
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
