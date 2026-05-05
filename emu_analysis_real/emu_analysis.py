"""
Parse EMU alignment outputs and generate annotated alignment tables.

Usage:
    python emu_analysis.py <emu_alignments.sam> <rel-abundance.tsv>

Arguments:
    emu_alignments.sam  - SAM file from EMU (requires --keep-files flag)
    rel-abundance.tsv   - Relative abundance output from EMU

Requires emu_taxonomy.tsv in the same directory as this script.

Outputs (written to current directory):
    emu_alignments_annotated.tsv - Full alignment table with species annotations
    emu_best_hit_per_read.tsv    - Best alignment per read
    emu_species_summary.tsv      - Per-species summary statistics
"""
import pandas as pd
import sys
import os
from collections import defaultdict
from sys import argv

if len(argv) < 3:
    print(__doc__)
    sys.exit(1)

emu_sam_file  = argv[1]
emu_rel_abund = argv[2]
emu_tax_file  = os.path.join(os.path.dirname(__file__), "..", "emu_taxonomy.tsv")


# ── Load taxonomy: tax_id -> species ─────────────────────────────────────────
tax_df  = pd.read_csv(emu_tax_file, sep="\t", dtype=str)
tax_map = dict(zip(tax_df["tax_id"], tax_df["species"]))

# ── Load EMU called species: set of tax_ids in the abundance output ───────────
abund_df     = pd.read_csv(emu_rel_abund, sep="\t", dtype=str)
called_taxids = set(abund_df["tax_id"].dropna())

# ── Parse SAM ────────────────────────────────────────────────────────────────
def parse_tags(fields):
    tags = {}
    for f in fields:
        if ":" in f:
            parts = f.split(":")
            if len(parts) == 3:
                tags[parts[0]] = parts[2]
    return tags

records = []
with open(emu_sam_file) as fh:
    for line in fh:
        if line.startswith("@"):
            continue
        cols = line.rstrip("\n").split("\t")
        read_id  = cols[0]
        flag     = int(cols[1])
        ref_name = cols[2]
        mapq     = int(cols[4])

        if ref_name == "*":   # unmapped
            continue

        # Parse reference: tax_id:emu_db:strain_id
        parts = ref_name.split(":")
        tax_id   = parts[0]
        strain_id = parts[2] if len(parts) >= 3 else None

        tags = parse_tags(cols[11:])
        AS   = int(tags.get("AS", -1))
        NM   = int(tags.get("NM", -1))
        de   = float(tags.get("de", float("nan")))
        is_primary = (flag & 256) == 0   # FLAG 256 = secondary

        species   = tax_map.get(tax_id, "unknown")
        in_emu    = tax_id in called_taxids

        records.append({
            "read_id":    read_id,
            "ref_name":   ref_name,
            "tax_id":     tax_id,
            "strain_id":  strain_id,
            "species":    species,
            "AS":         AS,
            "NM":         NM,
            "de":         de,
            "mapq":       mapq,
            "is_primary": is_primary,
            "in_emu":     in_emu,
        })

df = pd.DataFrame(records)

# ── Compute delta_AS: per-read difference from that read's best AS ────────────
best_as = df.groupby("read_id")["AS"].transform("max")
df["delta_AS"] = best_as - df["AS"]

# ── Write full alignment table ────────────────────────────────────────────────
out_all = "emu_alignments_annotated.tsv"
df.to_csv(out_all, sep="\t", index=False)
print(f"Written: {out_all}  ({len(df)} alignments, {df['read_id'].nunique()} reads)")

# ── Summary: best-hit per read (primary or max AS) ───────────────────────────
best = df.sort_values("AS", ascending=False).drop_duplicates("read_id").copy()

out_best = "emu_best_hit_per_read.tsv"
best.to_csv(out_best, sep="\t", index=False)
print(f"Written: {out_best}  ({len(best)} reads)")

# ── Distribution of divergence (de) for best hits ─────────────────────────────
print("\n--- Divergence (de) distribution for best hit per read ---")
print(best["de"].describe().round(4).to_string())

print("\n--- By in_emu status ---")
print(best.groupby("in_emu")["de"].describe().round(4).to_string())

# ── How many best-hit species end up called by EMU? ──────────────────────────
print("\n--- Best-hit species: called vs not called by EMU ---")
summary = (best.groupby(["species", "tax_id", "in_emu"])
               .agg(n_reads=("read_id", "count"),
                    mean_de=("de", "mean"),
                    mean_AS=("AS", "mean"))
               .reset_index()
               .sort_values("n_reads", ascending=False))
print(summary.to_string(index=False))

out_summary = "emu_species_summary.tsv"
summary.to_csv(out_summary, sep="\t", index=False)
print(f"\nWritten: {out_summary}")

# ── Delta AS distribution: how spread are alignments within a read? ───────────
print("\n--- delta_AS distribution (all alignments) ---")
print(df["delta_AS"].describe().round(1).to_string())
print(f"\nAlignments with delta_AS == 0 (tied-best): {(df['delta_AS'] == 0).sum()}")
