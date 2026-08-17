# Savont sensitivity to simulated read accuracy

This directory contains the read-accuracy sensitivity analysis. The design
matches the existing minimum-cluster-size and primary-clustering-threshold
panels: Savont is evaluated at 300, 1,318, 5,792, 25,450, and 111,828 input
reads, and a reported ASV is correct only when its best ground-truth alignment
has `NM=0`. Sensitivity is the fraction of all 27 Zymo V1--V8 ground-truth ASVs
recovered with an `NM=0` alignment. Terminal clipping is allowed, consistent
with the other parameter sweeps.

## Simulation design

- Badread 0.4.1 mean identities: 96%, 97%, 98%, and 99%.
- Maximum identity: 99.99%; identity standard deviation: 2.5 percentage points.
- Error and q-score models: `nanopore2023`.
- The requested fragment length is placed above the 1.36--1.41-kb linear
  reference lengths, causing Badread to simulate complete V1--V8 amplicons.
- Chimeras, random/junk reads, adapters, and glitches are disabled so that the
  sweep isolates read accuracy.
- Relative reference weights come from the mapping-derived mean-depth estimates
  in
  `results-FINAL-preprint/coverage/v1v8_zmock1_frac300000.0000000001.tsv`.
  Badread selects sequence proportional to `depth * reference length`, which is
  proportional to mapped read support. The exact weights are recorded in
  `ground_truth_abundance_weights.tsv`.
- Each accuracy has one deterministic maximum-depth simulation. The five input
  sizes are nested prefixes of that simulation, making depth comparisons within
  an accuracy level deterministic and directly analogous to nested subsampling.

The simulated depth is the number of raw input reads. Savont retains its default
read-quality filter; therefore, decreased sensitivity at lower simulated
accuracies can reflect both more error-disrupted sequence evidence and fewer
reads passing the default quality threshold. The empirical identity distribution
and fraction of reads with true identity at least 98% are reported for every
cell in `simulation_accuracy_audit.tsv`.

## Key outputs

- `read_accuracy_sensitivity.tsv`: precision and sensitivity for every grid cell.
- `read_accuracy_sensitivity_heatmaps.pdf` and `.svg`: manuscript-style panel.
- `simulation_accuracy_audit.tsv`: requested and empirical identity summaries.
- `simulation_realized_abundance.tsv`: realized versus mapping-derived abundance
  for every truth ASV and nested sample.
- `truth_asv_recovery.tsv`: exact recovery status for every truth ASV and grid
  cell.
- `ground_truth_abundance_weights.tsv`: input mapping-depth weights.
- `software_versions.txt`: exact tool versions.
- `logs/`, `results/`, and `alignments/`: complete run evidence.

## Reproduction

From the repository root:

```bash
MPLCONFIGDIR=/tmp/matplotlib-read-accuracy \
python reviewer_analyses/read_accuracy_sensitivity.py
```

The script reuses complete simulations and successful Savont runs, but stops for
inspection rather than silently overwriting a partial output.

## Results

All 20 Savont runs completed successfully. Precision was 100% in every one of
the 19 cells that produced at least one ASV; the 96%-accuracy, 300-read cell
reported no ASVs, so its precision is undefined rather than zero. Sensitivity
increased strongly with accuracy at lower read counts:

| Input reads | 96% | 97% | 98% | 99% |
|---:|---:|---:|---:|---:|
| 300 | 0.0% | 3.7% | 14.8% | 25.9% |
| 1,318 | 22.2% | 40.7% | 63.0% | 74.1% |
| 5,792 | 48.1% | 77.8% | 85.2% | 88.9% |
| 25,450 | 66.7% | 85.2% | 88.9% | 77.8% |
| 111,828 | 85.2% | 85.2% | 88.9% | 85.2% |

The empirical mean identities in the maximum-depth samples were 95.96%,
96.97%, 97.95%, and 98.98%, with empirical standard deviations of 2.47--2.52
percentage points. At maximum depth, every realized reference fraction was
within 0.17 percentage points of its mapping-derived target. At 111,828 input
reads, Savont retained 18.7%, 39.2%, 62.9%,
and 84.9% after its default quality filter as mean accuracy increased
from 96% to 99%. Thus, much of the accuracy effect at sparse depth is mediated
by effective post-filter read support.

