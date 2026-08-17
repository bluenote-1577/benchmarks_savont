#!/usr/bin/env python3
"""Run the BaNaNA OTU workflow on an already Barrnap-extracted 18S FASTQ."""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO


SAMPLE = "sample"


def count_fastq(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for _ in handle) // 4


def count_fasta(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    return sum(1 for _ in SeqIO.parse(path, "fasta"))


def abundance_sum(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                total += int(fields[1])
    return total


def command_text(command: list[str]) -> str:
    return " ".join(shlex.quote(value) for value in command)


def run(
    command: list[str],
    *,
    cwd: Path,
    env: dict[str, str],
    log: Path,
    stdout: Path | None = None,
) -> None:
    # Some BaNaNA stages launch parallel child processes. Capture each stage in
    # an isolated file first so their inherited output descriptors cannot
    # overwrite buffered text from an earlier stage in the cumulative log.
    with tempfile.NamedTemporaryFile(dir=log.parent, delete=False) as stage_log:
        stage_log_path = Path(stage_log.name)
        stage_log.write(f"\n$ {command_text(command)}\n".encode())
        stage_log.flush()
        try:
            if stdout is None:
                subprocess.run(
                    command,
                    cwd=cwd,
                    env=env,
                    stdout=stage_log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            else:
                stdout.parent.mkdir(parents=True, exist_ok=True)
                with stdout.open("wb") as output_handle:
                    subprocess.run(
                        command,
                        cwd=cwd,
                        env=env,
                        stdout=output_handle,
                        stderr=stage_log,
                        check=True,
                    )
        finally:
            stage_log.flush()
            with log.open("ab") as log_handle, stage_log_path.open("rb") as captured:
                shutil.copyfileobj(captured, log_handle)
            stage_log_path.unlink()


def write_summary(path: Path, values: dict[str, object]) -> None:
    with path.open("w") as handle:
        handle.write("metric\tvalue\n")
        for key, value in values.items():
            handle.write(f"{key}\t{value}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fastq", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--pr2-db", required=True, type=Path)
    parser.add_argument("--banana-dir", required=True, type=Path)
    parser.add_argument("--tool-prefix", required=True, type=Path)
    parser.add_argument("--min-length", required=True, type=int)
    parser.add_argument("--max-length", required=True, type=int)
    parser.add_argument("--min-mean-quality", default=90.0, type=float)
    parser.add_argument("--min-polishing-quality", default=20, type=int)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()

    input_fastq = args.input_fastq.resolve()
    output_dir = args.output_dir.resolve()
    pr2_db = args.pr2_db.resolve()
    banana_dir = args.banana_dir.resolve()
    tool_bin = (args.tool_prefix / "bin").resolve()
    summary_path = args.summary.resolve()
    log_path = args.log.resolve()

    samples = output_dir / "samples"
    if samples.exists():
        shutil.rmtree(samples)
    samples.mkdir(parents=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        "BaNaNA adapted entry point: input was already Barrnap-extracted; "
        "barrnap and extracting_rrna.py were intentionally skipped.\n"
    )

    env = os.environ.copy()
    env["PATH"] = f"{tool_bin}:{env.get('PATH', '')}"
    python = str(tool_bin / "python")
    filtlong = str(tool_bin / "filtlong")
    nanoplot = str(tool_bin / "NanoPlot")
    vsearch = str(tool_bin / "vsearch")
    racon = str(tool_bin / "racon")

    filtered_fastq = samples / f"filtlong_{SAMPLE}.fastq"
    rrna_fasta = samples / f"rrna_extracted_{SAMPLE}.fasta"
    nanoplot_dir = samples / f"nanoplot_{SAMPLE}"
    threshold_file = samples / f"clust_file_{SAMPLE}.txt"
    error_clusters = samples / f"clusters_error_{SAMPLE}"
    consensus = samples / f"consensus_{SAMPLE}.fasta"
    minimap_dir = samples / f"minimap_out_{SAMPLE}"
    combined_paf = samples / f"minimap_out_all_{SAMPLE}.paf"
    racon_fasta = samples / f"racon_{SAMPLE}.fasta"
    named_fasta = samples / f"racon_name_{SAMPLE}.fasta"
    merged = samples / "merged.fasta"
    nonchim_ref = samples / "nonchim_db.fasta"
    nonchim_denovo = samples / "nonchim_db_dn.fasta"
    final_clusters = samples / "clusters_final"
    pre_otus = samples / "pre_otus.fasta"
    final_dir = samples / "final"
    final_otus = final_dir / "otus.fasta"
    abundance = samples / f"abundance_{SAMPLE}.tsv"

    run(
        [
            filtlong,
            "--min_length",
            str(args.min_length),
            "--max_length",
            str(args.max_length),
            "--min_mean_q",
            str(args.min_mean_quality),
            str(input_fastq),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
        stdout=filtered_fastq,
    )

    SeqIO.convert(filtered_fastq, "fastq", rrna_fasta, "fasta")
    filtered_reads = count_fastq(filtered_fastq)
    if filtered_reads == 0:
        final_dir.mkdir(parents=True)
        final_otus.touch()
        abundance.touch()
        write_summary(
            summary_path,
            {
                "input_reads": count_fastq(input_fastq),
                "filtlong_reads": 0,
                "barrnap_skipped_already_extracted": True,
                "final_otus": 0,
                "otu_abundance_sum": 0,
            },
        )
        return

    run(
        [
            nanoplot,
            "--fastq",
            str(filtered_fastq.relative_to(output_dir)),
            "--tsv_stats",
            "-t",
            str(args.threads),
            "--info_in_report",
            "-o",
            str(nanoplot_dir.relative_to(output_dir)),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )
    run(
        [
            python,
            str(banana_dir / "scripts" / "calculate_clustering_threshold.py"),
            "-s",
            str((nanoplot_dir / "NanoStats.txt").relative_to(output_dir)),
            "-e",
            str(banana_dir / "files" / "P_error_table.tsv"),
            "-o",
            str(threshold_file.relative_to(output_dir)),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )
    threshold = threshold_file.read_text().strip()

    error_clusters.mkdir()
    run(
        [
            vsearch,
            "--cluster_fast",
            str(rrna_fasta.relative_to(output_dir)),
            "--id",
            threshold,
            "--threads",
            str(args.threads),
            "--clusters",
            f"samples/clusters_error_{SAMPLE}/cluster",
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )
    clusters_before_min4 = len(list(error_clusters.glob("cluster*")))

    alignments = error_clusters / "alignments"
    alignments.mkdir()
    run(
        [
            python,
            str(banana_dir / "scripts" / "mafft_consensus.py"),
            "-i",
            f"samples/clusters_error_{SAMPLE}/",
            "-a",
            f"samples/clusters_error_{SAMPLE}/alignments/",
            "-t",
            str(args.threads),
            "-o",
            str(consensus.relative_to(output_dir)),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )
    clusters_after_min4 = len(
        [path for path in error_clusters.glob("cluster*") if path.is_file()]
    )

    if count_fasta(consensus) == 0:
        final_dir.mkdir(parents=True)
        final_otus.touch()
        abundance.touch()
        write_summary(
            summary_path,
            {
                "input_reads": count_fastq(input_fastq),
                "filtlong_reads": filtered_reads,
                "barrnap_skipped_already_extracted": True,
                "clustering_identity": threshold,
                "error_clusters_before_min_size_4": clusters_before_min4,
                "error_clusters_after_min_size_4": clusters_after_min4,
                "final_otus": 0,
                "otu_abundance_sum": 0,
            },
        )
        return

    minimap_dir.mkdir()
    run(
        [
            python,
            str(banana_dir / "scripts" / "minimap.py"),
            "-c",
            str(consensus.relative_to(output_dir)),
            "-cl",
            f"samples/clusters_error_{SAMPLE}/",
            "-t",
            str(args.threads),
            "-o",
            f"samples/minimap_out_{SAMPLE}/",
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )
    with combined_paf.open("wb") as output_handle:
        for paf in sorted(minimap_dir.glob("*")):
            with paf.open("rb") as input_handle:
                shutil.copyfileobj(input_handle, output_handle)

    run(
        [
            racon,
            str(rrna_fasta.relative_to(output_dir)),
            "-q",
            str(args.min_polishing_quality),
            "-w",
            "500",
            "-t",
            str(args.threads),
            str(combined_paf.relative_to(output_dir)),
            str(consensus.relative_to(output_dir)),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
        stdout=racon_fasta,
    )
    run(
        [
            python,
            str(banana_dir / "scripts" / "add_sample_id.py"),
            "-i",
            str(racon_fasta.relative_to(output_dir)),
            "-sn",
            SAMPLE,
            "-o",
            str(named_fasta.relative_to(output_dir)),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )
    shutil.copyfile(named_fasta, merged)

    run(
        [
            vsearch,
            "--uchime_ref",
            str(merged.relative_to(output_dir)),
            "--db",
            str(pr2_db),
            "--threads",
            str(args.threads),
            "--nonchimeras",
            str(nonchim_ref.relative_to(output_dir)),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )
    run(
        [
            vsearch,
            "--uchime2_denovo",
            str(nonchim_ref.relative_to(output_dir)),
            "--threads",
            str(args.threads),
            "--nonchimeras",
            str(nonchim_denovo.relative_to(output_dir)),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )

    final_clusters.mkdir()
    run(
        [
            vsearch,
            "--cluster_fast",
            str(nonchim_denovo.relative_to(output_dir)),
            "--id",
            "0.99",
            "--threads",
            str(args.threads),
            "--clusters",
            "samples/clusters_final/cluster",
            "--centroids",
            str(pre_otus.relative_to(output_dir)),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )
    final_dir.mkdir()
    run(
        [
            python,
            str(banana_dir / "scripts" / "remove_Nseqs.py"),
            "-i",
            str(pre_otus.relative_to(output_dir)),
            "-o",
            str(final_otus.relative_to(output_dir)),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )
    run(
        [
            python,
            str(banana_dir / "scripts" / "abundance.py"),
            "-otu",
            str(final_otus.relative_to(output_dir)),
            "-fclu",
            "samples/clusters_final/",
            "-eclu",
            f"samples/clusters_error_{SAMPLE}/",
            "-sn",
            SAMPLE,
            "-o",
            str(abundance.relative_to(output_dir)),
        ],
        cwd=output_dir,
        env=env,
        log=log_path,
    )

    write_summary(
        summary_path,
        {
            "input_reads": count_fastq(input_fastq),
            "filtlong_reads": filtered_reads,
            "barrnap_skipped_already_extracted": True,
            "clustering_identity": threshold,
            "error_clusters_before_min_size_4": clusters_before_min4,
            "error_clusters_after_min_size_4": clusters_after_min4,
            "racon_consensuses": count_fasta(racon_fasta),
            "reference_nonchimeras": count_fasta(nonchim_ref),
            "denovo_nonchimeras": count_fasta(nonchim_denovo),
            "final_otus": count_fasta(final_otus),
            "otu_abundance_sum": abundance_sum(abundance),
        },
    )


if __name__ == "__main__":
    main()
