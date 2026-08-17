#!/usr/bin/env python3
"""Quantify full-depth ONT recovery of abundance-supported HiFi ASVs."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "hifi_supported_recovery"
MINIMAP2 = Path("/homes9/jshaw/bin/minimap2")
SAMPLES = ["AD1", "Soil1", "WWTP1", "Zfecal1"]
DEPTH_RE = re.compile(r"_depth_([0-9]+)(?:_|\s|$)")


def read_records(path: Path) -> dict[str, tuple[int, str]]:
    records: dict[str, tuple[int, str]] = {}
    name: str | None = None
    depth: int | None = None
    sequence_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None and depth is not None:
                    records[name] = (depth, "".join(sequence_parts).upper())
                name = line[1:].strip().split()[0]
                match = DEPTH_RE.search(line)
                if match is None:
                    raise ValueError(f"No depth in FASTA header: {line.rstrip()}")
                depth = int(match.group(1))
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())
    if name is not None and depth is not None:
        records[name] = (depth, "".join(sequence_parts).upper())
    return records


def canonical(sequence: str) -> str:
    reverse_complement = sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    return min(sequence, reverse_complement)


def parse_best_paf(path: Path) -> dict[str, tuple[int, str]]:
    best: dict[str, tuple[int, str]] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            nm = int(next(tag.rsplit(":", 1)[1] for tag in fields[12:] if tag.startswith("NM:i:")))
            query, target = fields[0], fields[5]
            if query not in best or nm < best[query][0]:
                best[query] = (nm, target)
    return best


def main() -> None:
    alignments = OUT / "alignments"
    alignments.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for sample in SAMPLES:
        hifi = BASE / "morten-zotu-data" / "pacbio_asvs" / f"{sample}_savont.fa"
        ont = BASE / "results_real" / "savont" / f"{sample}_ont_frac1" / "final_asvs.fasta"
        paf = alignments / f"{sample}_hifi_to_ont.paf"
        with paf.open("w") as output:
            subprocess.run(
                [str(MINIMAP2), "--cs", "--secondary=no", str(ont), str(hifi)],
                stdout=output,
                check=True,
            )

        hifi_records = read_records(hifi)
        ont_records = read_records(ont)
        hifi_depths = {name: record[0] for name, record in hifi_records.items()}
        ont_depths = {name: record[0] for name, record in ont_records.items()}
        ont_sequences = {canonical(record[1]) for record in ont_records.values()}
        best = parse_best_paf(paf)
        hifi_total = sum(hifi_depths.values())
        ont_max_relative_abundance = max(ont_depths.values()) / sum(ont_depths.values())
        cutoff = ont_max_relative_abundance / 10
        hifi_max_relative_abundance = max(hifi_depths.values()) / hifi_total

        for asv, depth in hifi_depths.items():
            nm, target = best.get(asv, (None, None))
            relative_abundance = depth / hifi_total
            rows.append(
                {
                    "sample": sample,
                    "hifi_asv": asv,
                    "hifi_depth": depth,
                    "hifi_relative_abundance": relative_abundance,
                    "hifi_max_relative_abundance": hifi_max_relative_abundance,
                    "ont_max_relative_abundance": ont_max_relative_abundance,
                    "one_tenth_ont_max_cutoff": cutoff,
                    "abundance_group": ">1/10 ONT maximum" if relative_abundance > cutoff else "<=1/10 ONT maximum",
                    "best_ont_asv": target,
                    "best_NM": nm,
                    "exactly_recovered": nm == 0,
                    "literal_full_length_match": canonical(hifi_records[asv][1]) in ont_sequences,
                }
            )

    observations = pd.DataFrame(rows)
    summaries = []
    for keys, group in observations.groupby(["sample", "abundance_group"], sort=False):
        sample, abundance_group = keys
        summaries.append(
            {
                "sample": sample,
                "abundance_group": abundance_group,
                "hifi_asvs": len(group),
                "exactly_recovered": int(group["exactly_recovered"].sum()),
                "exact_recovery_pct": 100 * group["exactly_recovered"].mean(),
            }
        )
    for abundance_group, group in observations.groupby("abundance_group", sort=False):
        summaries.append(
            {
                "sample": "Pooled",
                "abundance_group": abundance_group,
                "hifi_asvs": len(group),
                "exactly_recovered": int(group["exactly_recovered"].sum()),
                "exact_recovery_pct": 100 * group["exactly_recovered"].mean(),
            }
        )

    summary = pd.DataFrame(summaries)
    observations.to_csv(OUT / "hifi_supported_recovery_observations.tsv", sep="\t", index=False)
    summary.to_csv(
        OUT / "hifi_supported_recovery_summary.tsv", sep="\t", index=False, float_format="%.4f"
    )

    threshold_rows = []
    for sample, sample_data in observations.groupby("sample", sort=False):
        for basis, maximum_column in [
            ("ONT maximum relative abundance", "ont_max_relative_abundance"),
            ("HiFi maximum relative abundance", "hifi_max_relative_abundance"),
        ]:
            maximum = sample_data[maximum_column].iloc[0]
            for fraction in [1 / 3, 1 / 10]:
                for comparison, selected in [
                    (">", sample_data["hifi_relative_abundance"] > fraction * maximum),
                    ("<=", sample_data["hifi_relative_abundance"] <= fraction * maximum),
                ]:
                    group = sample_data[selected]
                    threshold_rows.append(
                        {
                            "sample": sample,
                            "threshold_basis": basis,
                            "threshold_fraction": fraction,
                            "comparison": comparison,
                            "hifi_asvs": len(group),
                            "NM0_recovered": int(group["exactly_recovered"].sum()),
                            "NM0_recovery_pct": 100 * group["exactly_recovered"].mean(),
                            "literal_full_length_matches": int(group["literal_full_length_match"].sum()),
                            "literal_full_length_match_pct": 100 * group["literal_full_length_match"].mean(),
                        }
                    )

    threshold_summary = pd.DataFrame(threshold_rows)
    pooled_rows = []
    for keys, group in threshold_summary.groupby(
        ["threshold_basis", "threshold_fraction", "comparison"], sort=False
    ):
        basis, fraction, comparison = keys
        total = int(group["hifi_asvs"].sum())
        nm0 = int(group["NM0_recovered"].sum())
        literal = int(group["literal_full_length_matches"].sum())
        pooled_rows.append(
            {
                "sample": "Pooled",
                "threshold_basis": basis,
                "threshold_fraction": fraction,
                "comparison": comparison,
                "hifi_asvs": total,
                "NM0_recovered": nm0,
                "NM0_recovery_pct": 100 * nm0 / total,
                "literal_full_length_matches": literal,
                "literal_full_length_match_pct": 100 * literal / total,
            }
        )
    threshold_summary = pd.concat([threshold_summary, pd.DataFrame(pooled_rows)], ignore_index=True)
    threshold_summary.to_csv(
        OUT / "hifi_supported_recovery_threshold_summary.tsv",
        sep="\t",
        index=False,
        float_format="%.6f",
    )
    print(summary.to_string(index=False))
    print("\nThreshold checks:")
    print(threshold_summary[threshold_summary["sample"] == "Pooled"].to_string(index=False))


if __name__ == "__main__":
    main()
