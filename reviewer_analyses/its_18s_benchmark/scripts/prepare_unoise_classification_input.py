#!/usr/bin/env python3
"""Wrap UNOISE zOTUs and OTU-table depths as a Savont classification input."""

from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path):
    header = None
    sequence: list[str] = []
    with path.open() as handle:
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


def read_counts(path: Path) -> tuple[str, dict[str, int]]:
    sample = "sample"
    counts: dict[str, int] = {}
    if not path.exists() or path.stat().st_size == 0:
        return sample, counts
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if fields[0].startswith("#OTU ID"):
                if len(fields) > 1:
                    sample = fields[1]
                continue
            if not fields or fields[0].startswith("#"):
                continue
            counts[fields[0]] = sum(int(float(value)) for value in fields[1:] if value)
    return sample, counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zotus", required=True, type=Path)
    parser.add_argument("--otutab", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample-name", required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    _, counts = read_counts(args.otutab)
    records = list(read_fasta(args.zotus)) if args.zotus.exists() else []
    renamed: list[tuple[str, str, int]] = []
    for zotu, sequence in records:
        depth = max(1, counts.get(zotu, 1))
        renamed.append((f"{zotu}_depth_{depth}", sequence, depth))

    with (args.output_dir / "final_asvs.fasta").open("w") as fasta:
        for header, sequence, _ in renamed:
            fasta.write(f">{header}\n{sequence}\n")
    with (args.output_dir / "feature-table.tsv").open("w") as table:
        table.write(f"#OTU ID\t{args.sample_name}\n")
        for header, _, depth in renamed:
            table.write(f"{header}\t{depth}\n")


if __name__ == "__main__":
    main()

