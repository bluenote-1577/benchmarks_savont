#!/usr/bin/env python3
"""Score expected-genus recovery, using binomial species names as aliases."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def normalize(value: str) -> str:
    return re.sub(r"\s+", " ", value.replace("_", " ").strip()).casefold()


def matches(alias: str, species: str, genus: str) -> bool:
    alias = normalize(alias)
    species = normalize(species)
    genus = normalize(genus)
    return genus == alias or species == alias or species.startswith(alias + " ")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-table", required=True, type=Path)
    parser.add_argument("--expected", required=True, type=Path)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--reads", required=True, type=int)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--details", required=True, type=Path)
    args = parser.parse_args()

    expected = []
    with args.expected.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["marker"] == args.marker:
                expected.append((row["canonical_genus"], row["aliases"].split(";")))

    classified = []
    with args.species_table.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            classified.append(row)

    detail_rows = []
    for canonical, aliases in expected:
        hits = [
            row
            for row in classified
            if any(matches(alias, row.get("species", ""), row.get("genus", "")) for alias in aliases)
        ]
        detail_rows.append(
            {
                "marker": args.marker,
                "method": args.method,
                "reads": args.reads,
                "expected_genus": canonical,
                "recovered": bool(hits),
                "matched_species": "|".join(sorted({row.get("species", "") for row in hits})),
                "matched_genera": "|".join(sorted({row.get("genus", "") for row in hits})),
                "summed_abundance": sum(float(row.get("abundance", 0)) for row in hits),
            }
        )

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w") as handle:
        columns = ["marker", "method", "reads", "genera_recovered", "expected_genera"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "marker": args.marker,
                "method": args.method,
                "reads": args.reads,
                "genera_recovered": sum(row["recovered"] for row in detail_rows),
                "expected_genera": len(detail_rows),
            }
        )
    with args.details.open("w") as handle:
        columns = list(detail_rows[0])
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(detail_rows)


if __name__ == "__main__":
    main()

