#!/usr/bin/env python3
"""Slice the best nhmmer hit per read from FASTQ while retaining qualities."""

from __future__ import annotations

import argparse
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def load_hits(path: Path) -> dict[str, tuple[int, int, str, float]]:
    hits: dict[str, tuple[int, int, str, float]] = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            target = fields[0]
            alignment_from, alignment_to = int(fields[6]), int(fields[7])
            strand = fields[11]
            evalue = float(fields[12])
            start, end = sorted((alignment_from, alignment_to))
            previous = hits.get(target)
            if previous is None or evalue < previous[3]:
                hits[target] = (start, end, strand, evalue)
    return hits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tblout", required=True, type=Path)
    parser.add_argument("--fastq", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    hits = load_hits(args.tblout)
    written = 0
    with args.fastq.open() as source, args.output.open("w") as output:
        while True:
            header = source.readline()
            if not header:
                break
            sequence = source.readline().rstrip("\r\n")
            plus = source.readline()
            quality = source.readline().rstrip("\r\n")
            if not plus or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record after {header.rstrip()}")
            read_id = header[1:].split()[0]
            hit = hits.get(read_id)
            if hit is None:
                continue
            start, end, strand, _ = hit
            subsequence = sequence[start - 1 : end]
            subquality = quality[start - 1 : end]
            description = f"18S_rRNA:{start}-{end}"
            if strand == "-":
                subsequence = reverse_complement(subsequence)
                subquality = subquality[::-1]
                description += " rc"
            output.write(
                f"@{read_id} {description}\n{subsequence}\n+\n{subquality}\n"
            )
            written += 1
    print(f"{written} / {len(hits)} hit reads written to {args.output}")


if __name__ == "__main__":
    main()
