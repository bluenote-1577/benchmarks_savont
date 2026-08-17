#!/usr/bin/env python3
"""Combine per-run taxon-recovery tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def combine(paths: list[Path], output: Path) -> None:
    tables = [pd.read_csv(path, sep="\t") for path in paths]
    combined = pd.concat(tables, ignore_index=True).sort_values(["marker", "method", "reads"])
    output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output, sep="\t", index=False, float_format="%.6f")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summaries", nargs="+", required=True, type=Path)
    parser.add_argument("--details", nargs="+", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    parser.add_argument("--details-output", required=True, type=Path)
    args = parser.parse_args()
    combine(args.summaries, args.summary_output)
    combine(args.details, args.details_output)


if __name__ == "__main__":
    main()

