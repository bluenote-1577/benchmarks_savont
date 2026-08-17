#!/usr/bin/env python3
"""Convert BaNaNA OTUs and abundances into a Savont classification input."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO


def read_abundances(path: Path) -> dict[str, int]:
    abundances: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                abundances[fields[0]] = int(fields[1])
    return abundances


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--otus", required=True, type=Path)
    parser.add_argument("--abundance", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample-name", required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    abundances = read_abundances(args.abundance)
    records = list(SeqIO.parse(args.otus, "fasta")) if args.otus.stat().st_size else []

    missing = [record.id for record in records if abundances.get(record.id, 0) <= 0]
    if missing:
        raise ValueError(f"BaNaNA OTUs lack positive abundance estimates: {missing}")

    fasta_path = args.output_dir / "final_asvs.fasta"
    table_path = args.output_dir / "feature-table.tsv"
    map_path = args.output_dir / "otu_id_map.tsv"
    with fasta_path.open("w") as fasta, table_path.open("w") as table, map_path.open("w") as id_map:
        table.write(f"#OTU ID\t{args.sample_name}\n")
        id_map.write("banana_otu\tsavont_header\tdepth\n")
        for index, record in enumerate(records, start=1):
            depth = abundances[record.id]
            header = f"BaNaNA_OTU{index}_depth_{depth}"
            fasta.write(f">{header}\n{str(record.seq)}\n")
            table.write(f"{header}\t{depth}\n")
            id_map.write(f"{record.id}\t{header}\t{depth}\n")


if __name__ == "__main__":
    main()
