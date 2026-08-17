#!/usr/bin/env python3
"""Create uniform, deterministic, nested FASTQ subsamples in one input pass."""

from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def percentile(sorted_values: list[int], fraction: float) -> float:
    if not sorted_values:
        return float("nan")
    position = fraction * (len(sorted_values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    weight = position - lower
    return sorted_values[lower] * (1 - weight) + sorted_values[upper] * weight


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--depths", required=True, nargs="+", type=int)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--min-length", required=True, type=int)
    parser.add_argument("--max-length", required=True, type=int)
    parser.add_argument("--audit", required=True, type=Path)
    args = parser.parse_args()

    depths = sorted(set(args.depths))
    maximum = max(depths)
    rng = random.Random(args.seed)
    reservoir: list[tuple[str, str, str, str]] = []
    lengths: list[int] = []
    raw_records = 0
    nonempty_records = 0
    zero_length_records = 0
    length_eligible_records = 0

    with open_text(args.input) as handle:
        while True:
            header = handle.readline()
            if header == "":
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            raw_records += 1
            if sequence == "" or plus == "" or quality == "":
                raise ValueError(f"Truncated FASTQ record {raw_records} in {args.input}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record {raw_records} in {args.input}")
            sequence_length = len(sequence.rstrip("\r\n"))
            quality_length = len(quality.rstrip("\r\n"))
            if sequence_length != quality_length:
                raise ValueError(f"Sequence/quality length mismatch at record {raw_records}")
            lengths.append(sequence_length)
            if sequence_length == 0:
                zero_length_records += 1
                continue
            nonempty_records += 1
            if args.min_length <= sequence_length <= args.max_length:
                length_eligible_records += 1

            record = (header, sequence, plus, quality)
            if len(reservoir) < maximum:
                reservoir.append(record)
            else:
                replacement = rng.randrange(nonempty_records)
                if replacement < maximum:
                    reservoir[replacement] = record

    if nonempty_records < maximum:
        raise ValueError(
            f"Requested {maximum} reads but only {nonempty_records} nonempty records are available"
        )
    rng.shuffle(reservoir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for depth in depths:
        output = args.output_dir / f"{args.prefix}_{depth}.fastq"
        with output.open("w") as handle:
            for record in reservoir[:depth]:
                handle.writelines(record)

    sorted_lengths = sorted(lengths)
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    with args.audit.open("w") as handle:
        handle.write(
            "marker\tsource\traw_records\tnonempty_records\tzero_length_records\t"
            "savont_length_eligible\tconfigured_min_length\tconfigured_max_length\t"
            "observed_min_length\tq01_length\tq05_length\tq25_length\tmedian_length\t"
            "q75_length\tq95_length\tq99_length\tobserved_max_length\tseed\n"
        )
        values = [
            args.prefix,
            str(args.input),
            raw_records,
            nonempty_records,
            zero_length_records,
            length_eligible_records,
            args.min_length,
            args.max_length,
            min(sorted_lengths),
            percentile(sorted_lengths, 0.01),
            percentile(sorted_lengths, 0.05),
            percentile(sorted_lengths, 0.25),
            percentile(sorted_lengths, 0.50),
            percentile(sorted_lengths, 0.75),
            percentile(sorted_lengths, 0.95),
            percentile(sorted_lengths, 0.99),
            max(sorted_lengths),
            args.seed,
        ]
        handle.write("\t".join(map(str, values)) + "\n")


if __name__ == "__main__":
    main()
