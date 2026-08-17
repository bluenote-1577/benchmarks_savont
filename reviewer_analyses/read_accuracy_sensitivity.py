#!/usr/bin/env python3
"""Benchmark Savont across Badread-simulated ONT accuracies and read depths.

The design mirrors the two existing referee parameter sweeps: five nested read
depths, exact (NM=0) matching to all 27 Zymo V1-V8 ground-truth ASVs, and two
precision/sensitivity heatmaps.  Ground-truth abundance weights come from the
300,000-read mapping-derived coverage table used by the main benchmark.
"""

from __future__ import annotations

import csv
import gzip
import re
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "read_accuracy_sensitivity"
SAVONT = Path("/homes9/jshaw/.cargo/bin/savont")
MINIMAP2 = Path("/homes9/jshaw/bin/minimap2")
BADREAD = Path("/homes9/jshaw/miniforge3/bin/badread")
TRUTH = BASE / "morten_data" / "V1V8_bacterial_amplicons.fasta.gz"
COVERAGE = (
    BASE
    / "results-FINAL-preprint"
    / "coverage"
    / "v1v8_zmock1_frac300000.0000000001.tsv"
)
MEAN_ACCURACIES = [96, 97, 98, 99]
MAX_ACCURACY = 99.99
ACCURACY_SD = 2.5
READ_DEPTHS = [300, 1318, 5792, 25450, 111828]
SIMULATED_BASES = "160M"
SEED = 727271
THREADS_PER_RUN = 10
PARALLEL_RUNS = 2
IDENTITY_RE = re.compile(r"read_identity=([0-9.]+)%")
READ_FILTER_RE = re.compile(
    r"Number of valid reads\s+-\s+(\d+)\. Number of reads below quality threshold - (\d+)\."
)


def fasta_records(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header: str | None = None
    sequence: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:].split()[0]
                sequence = []
            elif header is not None:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def count_fasta(path: Path) -> int:
    return sum(1 for _ in fasta_records(path))


def count_fastq(path: Path) -> int:
    with path.open() as handle:
        lines = sum(1 for _ in handle)
    if lines % 4:
        raise ValueError(f"FASTQ does not contain complete four-line records: {path}")
    return lines // 4


def build_weighted_reference() -> dict[str, float]:
    OUT.mkdir(parents=True, exist_ok=True)
    weighted = OUT / "ground_truth_weighted.fasta"
    audit = OUT / "ground_truth_abundance_weights.tsv"
    with COVERAGE.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    coverage_rows = {row["rname"]: row for row in rows}
    records = list(fasta_records(TRUTH))
    truth_names = {name for name, _ in records}
    if truth_names != set(coverage_rows):
        raise ValueError("Coverage-table targets do not exactly match the truth FASTA")

    weighted_support = {
        name: float(coverage_rows[name]["meandepth"]) * len(sequence)
        for name, sequence in records
    }
    total_support = sum(weighted_support.values())
    expected_fractions = {name: value / total_support for name, value in weighted_support.items()}

    with weighted.open("w") as output:
        for name, sequence in records:
            depth = float(coverage_rows[name]["meandepth"])
            output.write(f">{name} depth={depth:.6f}\n{sequence}\n")
    with audit.open("w") as output:
        output.write(
            "truth_asv\tlength\tmapped_reads\tmean_mapped_depth\t"
            "badread_selection_weight\texpected_read_fraction\n"
        )
        for name, sequence in records:
            row = coverage_rows[name]
            output.write(
                f"{name}\t{len(sequence)}\t{row['numreads']}\t{float(row['meandepth']):.6f}\t"
                f"{weighted_support[name]:.6f}\t{expected_fractions[name]:.9f}\n"
            )
    return expected_fractions


def simulate_one(mean_accuracy: int) -> Path:
    simulated_dir = OUT / "simulated"
    simulated_dir.mkdir(parents=True, exist_ok=True)
    output = simulated_dir / f"mean_accuracy_{mean_accuracy}.fastq"
    log = OUT / "logs" / f"badread_mean_accuracy_{mean_accuracy}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and count_fastq(output) >= max(READ_DEPTHS):
        return output
    if output.exists() and output.stat().st_size:
        raise RuntimeError(f"Incomplete simulation requires inspection: {output}")

    command = [
        str(BADREAD),
        "simulate",
        "--reference",
        str(OUT / "ground_truth_weighted.fasta"),
        "--quantity",
        SIMULATED_BASES,
        "--length",
        "2000,1",
        "--identity",
        f"{mean_accuracy},{MAX_ACCURACY},{ACCURACY_SD}",
        "--error_model",
        "nanopore2023",
        "--qscore_model",
        "nanopore2023",
        "--seed",
        str(SEED + mean_accuracy * 1000),
        "--start_adapter",
        "0,0",
        "--end_adapter",
        "0,0",
        "--junk_reads",
        "0",
        "--random_reads",
        "0",
        "--chimeras",
        "0",
        "--glitches",
        "0,0,0",
    ]
    with output.open("w") as fastq, log.open("w") as log_handle:
        subprocess.run(command, stdout=fastq, stderr=log_handle, check=True)
    observed = count_fastq(output)
    if observed < max(READ_DEPTHS):
        raise RuntimeError(f"Badread produced only {observed:,} reads for {mean_accuracy}%")
    return output


def make_nested_subsamples(
    mean_accuracy: int, expected_fractions: dict[str, float]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    source = OUT / "simulated" / f"mean_accuracy_{mean_accuracy}.fastq"
    subsample_dir = OUT / "subsampled"
    subsample_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        depth: subsample_dir / f"mean_accuracy_{mean_accuracy}_reads_{depth}.fastq"
        for depth in READ_DEPTHS
    }
    audit_rows: list[dict[str, object]] = []
    abundance_rows: list[dict[str, object]] = []

    handles = {depth: path.open("w") for depth, path in paths.items()}
    identities = {depth: [] for depth in READ_DEPTHS}
    taxa = {depth: Counter() for depth in READ_DEPTHS}
    try:
        with source.open() as input_handle:
            for index in range(1, max(READ_DEPTHS) + 1):
                record = [input_handle.readline() for _ in range(4)]
                if any(line == "" for line in record):
                    raise RuntimeError(f"Simulation ended before {max(READ_DEPTHS):,} reads: {source}")
                header = record[0].rstrip("\r\n")
                match = IDENTITY_RE.search(header)
                if match is None:
                    raise ValueError(f"Badread header lacks read_identity: {header}")
                identity = float(match.group(1))
                fields = header.split()
                if len(fields) < 2:
                    raise ValueError(f"Badread header lacks a source target: {header}")
                source_asv = fields[1].split(",", 1)[0]
                for depth in READ_DEPTHS:
                    if index <= depth:
                        handles[depth].writelines(record)
                        identities[depth].append(identity)
                        taxa[depth][source_asv] += 1
    finally:
        for handle in handles.values():
            handle.close()

    for depth in READ_DEPTHS:
        values = np.asarray(identities[depth], dtype=float)
        audit_rows.append(
            {
                "requested_mean_accuracy_pct": mean_accuracy,
                "maximum_accuracy_pct": MAX_ACCURACY,
                "requested_accuracy_sd_pct": ACCURACY_SD,
                "read_depth": depth,
                "actual_reads": len(values),
                "empirical_mean_identity_pct": values.mean(),
                "empirical_sd_identity_pct": values.std(ddof=1) if len(values) > 1 else np.nan,
                "empirical_min_identity_pct": values.min(),
                "empirical_max_identity_pct": values.max(),
                "fraction_reads_identity_ge_98": np.mean(values >= 98.0),
            }
        )
        for truth_asv in expected_fractions:
            count = taxa[depth][truth_asv]
            abundance_rows.append(
                {
                    "requested_mean_accuracy_pct": mean_accuracy,
                    "read_depth": depth,
                    "truth_asv": truth_asv,
                    "simulated_reads": count,
                    "realized_read_fraction": count / depth,
                    "mapping_derived_expected_fraction": expected_fractions[truth_asv],
                }
            )
    return audit_rows, abundance_rows


def is_complete(run_dir: Path) -> bool:
    fasta = run_dir / "final_asvs.fasta"
    if not fasta.exists():
        return False
    logs = sorted(run_dir.glob("savont_*.log"))
    return bool(logs) and "COMPLETED SUCCESSFULLY" in logs[-1].read_text(errors="replace")


def run_savont(mean_accuracy: int, read_depth: int) -> tuple[int, int]:
    reads = OUT / "subsampled" / f"mean_accuracy_{mean_accuracy}_reads_{read_depth}.fastq"
    run_dir = OUT / "results" / f"mean_accuracy_{mean_accuracy}" / f"reads_{read_depth}"
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    if is_complete(run_dir):
        return mean_accuracy, read_depth
    if run_dir.exists() and any(run_dir.iterdir()):
        raise RuntimeError(f"Incomplete non-empty run directory requires inspection: {run_dir}")
    subprocess.run(
        [
            str(SAVONT),
            "asv",
            "-t",
            str(THREADS_PER_RUN),
            str(reads),
            "-o",
            str(run_dir),
        ],
        check=True,
    )
    if not is_complete(run_dir):
        raise RuntimeError(f"Savont exited without a successful completion marker: {run_dir}")
    return mean_accuracy, read_depth


def parse_nm(fields: list[str]) -> int:
    for tag in fields[12:]:
        if tag.startswith("NM:i:"):
            return int(tag.rsplit(":", 1)[1])
    raise ValueError("PAF record is missing NM:i")


def score_one(mean_accuracy: int, read_depth: int, truth_count: int) -> dict[str, object]:
    run_dir = OUT / "results" / f"mean_accuracy_{mean_accuracy}" / f"reads_{read_depth}"
    asvs = run_dir / "final_asvs.fasta"
    paf = OUT / "alignments" / f"mean_accuracy_{mean_accuracy}_reads_{read_depth}.paf"
    paf.parent.mkdir(parents=True, exist_ok=True)
    with paf.open("w") as output:
        subprocess.run(
            [str(MINIMAP2), "--cs", "--secondary=no", str(TRUTH), str(asvs)],
            stdout=output,
            check=True,
        )

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
    precision = exact_asv_count / asv_count if asv_count else np.nan
    sensitivity = len(exact_targets) / truth_count
    logs = sorted(run_dir.glob("savont_*.log"))
    if not logs:
        raise RuntimeError(f"Missing Savont log: {run_dir}")
    filter_match = READ_FILTER_RE.search(logs[-1].read_text(errors="replace"))
    if filter_match is None:
        raise ValueError(f"Could not parse Savont read-filter counts: {logs[-1]}")
    length_valid_reads, reads_below_quality = map(int, filter_match.groups())
    return {
        "mean_read_accuracy_pct": mean_accuracy,
        "read_depth": read_depth,
        "length_valid_reads": length_valid_reads,
        "reads_below_quality_threshold": reads_below_quality,
        "reads_retained_after_quality": length_valid_reads - reads_below_quality,
        "fraction_reads_retained_after_quality": (
            (length_valid_reads - reads_below_quality) / length_valid_reads
            if length_valid_reads
            else np.nan
        ),
        "asv_count": asv_count,
        "exact_asv_count": exact_asv_count,
        "false_positive_count": asv_count - exact_asv_count,
        "truth_asvs_recovered": len(exact_targets),
        "truth_asv_count": truth_count,
        "precision": precision,
        "sensitivity": sensitivity,
    }


def truth_recovery_rows(mean_accuracy: int, read_depth: int) -> list[dict[str, object]]:
    paf = OUT / "alignments" / f"mean_accuracy_{mean_accuracy}_reads_{read_depth}.paf"
    exact_targets: set[str] = set()
    with paf.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if parse_nm(fields) == 0:
                exact_targets.add(fields[5])
    return [
        {
            "mean_read_accuracy_pct": mean_accuracy,
            "read_depth": read_depth,
            "truth_asv": name,
            "recovered_nm0": name in exact_targets,
        }
        for name, _ in fasta_records(TRUTH)
    ]


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
    for ax, metric, title in zip(
        axes, ["precision", "sensitivity"], ["Precision", "Sensitivity"]
    ):
        matrix = (
            results.pivot(index="read_depth", columns="mean_read_accuracy_pct", values=metric)
            .reindex(index=READ_DEPTHS, columns=MEAN_ACCURACIES)
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
        ax.set_xlabel("Mean simulated read accuracy (%)")
        ax.set_xticks(range(len(MEAN_ACCURACIES)), MEAN_ACCURACIES)
        ax.set_yticks(range(len(READ_DEPTHS)), [f"{depth:,}" for depth in READ_DEPTHS])
        ax.set_ylabel("Simulated ONT reads")
        for spine in ax.spines.values():
            spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=axes, fraction=0.025, pad=0.02)
    colorbar.set_ticks(np.arange(0, 1.01, 0.2))
    fig.savefig(OUT / "read_accuracy_sensitivity_heatmaps.svg", bbox_inches="tight")
    fig.savefig(OUT / "read_accuracy_sensitivity_heatmaps.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    expected_fractions = build_weighted_reference()
    versions = []
    for tool in (SAVONT, BADREAD, MINIMAP2):
        result = subprocess.run([str(tool), "--version"], capture_output=True, text=True, check=True)
        versions.append(f"{tool}: {(result.stdout or result.stderr).strip()}")
    (OUT / "software_versions.txt").write_text("\n".join(versions) + "\n")

    with ThreadPoolExecutor(max_workers=len(MEAN_ACCURACIES)) as executor:
        futures = {executor.submit(simulate_one, accuracy): accuracy for accuracy in MEAN_ACCURACIES}
        for future in as_completed(futures):
            path = future.result()
            print(f"Completed/read existing simulation: {path}", flush=True)

    audit_rows: list[dict[str, object]] = []
    abundance_rows: list[dict[str, object]] = []
    for accuracy in MEAN_ACCURACIES:
        audit, abundance = make_nested_subsamples(accuracy, expected_fractions)
        audit_rows.extend(audit)
        abundance_rows.extend(abundance)
    pd.DataFrame(audit_rows).to_csv(
        OUT / "simulation_accuracy_audit.tsv", sep="\t", index=False, float_format="%.6f"
    )
    pd.DataFrame(abundance_rows).to_csv(
        OUT / "simulation_realized_abundance.tsv", sep="\t", index=False, float_format="%.9f"
    )

    jobs = [(accuracy, depth) for accuracy in MEAN_ACCURACIES for depth in READ_DEPTHS]
    with ThreadPoolExecutor(max_workers=PARALLEL_RUNS) as executor:
        futures = {executor.submit(run_savont, *job): job for job in jobs}
        for future in as_completed(futures):
            accuracy, depth = future.result()
            print(f"Completed/read existing: accuracy={accuracy}%, reads={depth:,}", flush=True)

    truth_count = count_fasta(TRUTH)
    rows = [score_one(accuracy, depth, truth_count) for accuracy, depth in jobs]
    results = pd.DataFrame(rows).sort_values(["read_depth", "mean_read_accuracy_pct"])
    results.to_csv(
        OUT / "read_accuracy_sensitivity.tsv", sep="\t", index=False, float_format="%.6f"
    )
    recovery_rows = [
        row
        for accuracy, depth in jobs
        for row in truth_recovery_rows(accuracy, depth)
    ]
    pd.DataFrame(recovery_rows).to_csv(
        OUT / "truth_asv_recovery.tsv", sep="\t", index=False
    )
    plot_heatmaps(results)
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
