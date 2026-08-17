#!/usr/bin/env python3
"""Run UNOISE3, treating an empty dereplicated FASTA as zero recovered zOTUs."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--usearch", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--minsize", type=int)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    if not args.input.exists() or args.input.stat().st_size == 0:
        args.output.write_text("")
        args.log.write_text(
            "No dereplicated sequences met -minuniquesize 2; "
            "recording zero recovered zOTUs.\n"
        )
        return

    command = [str(args.usearch), "-unoise3", str(args.input), "-zotus", str(args.output)]
    if args.minsize is not None:
        command.extend(["-minsize", str(args.minsize)])
    with args.log.open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)


if __name__ == "__main__":
    main()
