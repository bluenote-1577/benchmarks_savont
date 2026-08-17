#!/usr/bin/env python3
"""Run Savont classification, with valid empty outputs for zero-ASV runs."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


SPECIES_HEADER = "abundance\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom\n"
GENUS_HEADER = "abundance\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom\n"
MAPPING_HEADER = (
    "asv_header\tdepth\talignment_identity\tnumber_mismatches\ttax_id\tspecies\tgenus\t"
    "family\torder\tclass\tphylum\tclade\tsuperkingdom\treference\n"
)


def has_sequences(path: Path) -> bool:
    with path.open() as handle:
        return any(line.startswith(">") for line in handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--savont", required=True, type=Path)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    fasta = args.input_dir / "final_asvs.fasta"

    if has_sequences(fasta):
        command = [
            str(args.savont),
            "classify",
            "-i",
            str(args.input_dir),
            "-o",
            str(args.output_dir),
            "-d",
            str(args.database),
            "-t",
            str(args.threads),
        ]
        with args.log.open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    else:
        (args.output_dir / "species_abundance.tsv").write_text(SPECIES_HEADER)
        (args.output_dir / "genus_abundance.tsv").write_text(GENUS_HEADER)
        (args.output_dir / "asv_mappings.tsv").write_text(MAPPING_HEADER)
        args.log.write_text("No ASVs were produced; wrote empty classification tables.\n")


if __name__ == "__main__":
    main()

