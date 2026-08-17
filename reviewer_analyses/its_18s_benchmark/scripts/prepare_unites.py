#!/usr/bin/env python3
"""Prepare Savont's short-header UNITE 10.0 database from the local archive."""

from __future__ import annotations

import argparse
import io
import tarfile
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fasta_output = args.output_dir / "unites_sequences.fasta"
    taxonomy_output = args.output_dir / "unites_taxonomy.tsv"

    with tarfile.open(args.archive, "r:gz") as archive:
        members = [
            member
            for member in archive.getmembers()
            if member.name.endswith(".fasta") and not member.name.endswith("_dev.fasta")
        ]
        if len(members) != 1:
            raise ValueError(f"Expected one primary FASTA in {args.archive}; found {len(members)}")
        binary = archive.extractfile(members[0])
        if binary is None:
            raise ValueError(f"Could not read {members[0].name} from {args.archive}")
        source = io.TextIOWrapper(binary)
        with fasta_output.open("w") as fasta, taxonomy_output.open("w") as taxonomy:
            taxonomy.write("sh_id\ttaxonomy\n")
            for line_number, line in enumerate(source, start=1):
                if not line.startswith(">"):
                    fasta.write(line)
                    continue
                header = line[1:].rstrip("\r\n")
                parts = header.split("|", 4)
                if len(parts) != 5:
                    raise ValueError(f"Malformed UNITE header at line {line_number}: {header[:120]}")
                sh_id, lineage = parts[2], parts[4]
                fasta.write(f">{sh_id}\n")
                taxonomy.write(f"{sh_id}\t{lineage}\n")

    (args.output_dir / ".savont_db").write_text("unites-10.0\n")


if __name__ == "__main__":
    main()

