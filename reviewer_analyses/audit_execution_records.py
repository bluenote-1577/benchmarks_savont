#!/usr/bin/env python3
"""Verify that the compact execution records expected in the Git copy exist.

This audits retained provenance records rather than rerunning the workflows.
Raw reads, databases, alignments, and full result directories remain outside
Git.
"""

from __future__ import annotations

import gzip
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def read_text(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def audit_logs(
    label: str,
    pattern: str,
    expected: int,
    required_text: str | None = None,
    allow_empty: bool = False,
) -> tuple[str, int, int, str]:
    paths = sorted(ROOT.glob(pattern))
    valid = 0
    for path in paths:
        if path.stat().st_size == 0 and not allow_empty:
            continue
        if required_text is None or required_text in read_text(path):
            valid += 1
    status = "PASS" if len(paths) == expected and valid == expected else "FAIL"
    return label, len(paths), valid, status


def main() -> int:
    checks = [
        audit_logs(
            "ITS/18S Savont runs",
            "its_18s_benchmark/logs/savont/*.log",
            12,
            "SAVONT COMPLETED SUCCESSFULLY",
        ),
        audit_logs(
            "ITS/18S UNOISE records",
            "its_18s_benchmark/logs/unoise/**/*.log",
            72,
            allow_empty=True,
        ),
        audit_logs(
            "ITS/18S classification records",
            "its_18s_benchmark/logs/classification/**/*.log",
            42,
        ),
        audit_logs(
            "18S BaNaNA command logs",
            "its_18s_benchmark/logs/banana/*.log.gz",
            6,
            "abundance.py",
        ),
        audit_logs(
            "18S BaNaNA stage summaries",
            "its_18s_benchmark/results/banana/18s_*/run_summary.tsv",
            6,
        ),
        audit_logs(
            "nested-subsampling audits",
            "its_18s_benchmark/results/subsampled/*.audit.tsv",
            2,
        ),
        audit_logs(
            "minimum-cluster Savont runs",
            "min_cluster_sensitivity/run_logs/**/*.log",
            45,
            "SAVONT COMPLETED SUCCESSFULLY",
        ),
        audit_logs(
            "primary-threshold Savont runs",
            "primary_clustering_threshold_sensitivity/run_logs/**/*.log",
            25,
            "SAVONT COMPLETED SUCCESSFULLY",
        ),
        audit_logs(
            "read-accuracy Savont runs",
            "read_accuracy_sensitivity/run_logs/savont/**/*.log",
            20,
            "SAVONT COMPLETED SUCCESSFULLY",
        ),
        audit_logs(
            "read-accuracy Badread simulations",
            "read_accuracy_sensitivity/run_logs/badread/*.log.gz",
            4,
            "Badread v0.4.1",
        ),
        audit_logs(
            "HiFi zOTU direct-mapping records",
            "zotu_recall_by_depth/logs/*.log",
            48,
        ),
    ]

    print("record_set\tfound\tvalid\tstatus")
    for row in checks:
        print("\t".join(map(str, row)))
    return int(any(row[-1] != "PASS" for row in checks))


if __name__ == "__main__":
    raise SystemExit(main())
