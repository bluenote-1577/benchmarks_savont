# Reviewer analyses

This directory contains result tables, plotting code, methods, and run records
for the manuscript revision.

## Regenerate all reviewer figures

Create the small plotting environment once:

```bash
conda env create -f reviewer_analyses/plotting_environment.yaml
conda activate savont-reviewer-plots
```

Then run from the repository root:

```bash
MPLCONFIGDIR=/tmp/matplotlib-reviewer \
python reviewer_analyses/plot_reviewer_figures.py
```

Use `--only` with one or more of `abundance`, `min-cluster`,
`primary-threshold`, `accuracy`, `zotu`, or `its-18s` to regenerate selected
figures.

## Plot inputs

| Figure | Compact input table(s) |
|---|---|
| HiFi agreement by Savont ASV abundance | `../real_notebook/asv_fidelity_by_abundance_subsamples.tsv` |
| Minimum cluster-size sensitivity | `min_cluster_sensitivity/min_cluster_sensitivity.tsv` |
| Primary-clustering-threshold sensitivity | `primary_clustering_threshold_sensitivity/primary_clustering_threshold_sensitivity.tsv` |
| Simulated read-accuracy sensitivity | `read_accuracy_sensitivity/read_accuracy_sensitivity.tsv` |
| HiFi zOTU recall by observed ONT depth | `zotu_recall_by_depth/zotu_recall_by_dataset.tsv` and `zotu_recall_by_depth/zotu_recall_by_depth_summary.tsv` |
| ITS/18S expected-genus recovery | `its_18s_benchmark/expected_genera_recovery.tsv` |

Additional observation, audit, and provenance tables support the plotted
summaries. The analysis scripts record how these tables were produced.
