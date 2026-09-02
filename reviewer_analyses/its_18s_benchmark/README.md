# ITS and 18S expected-genus recovery benchmark (supplementary)

This directory contains a reproducible comparison of Savont and UNOISE3 on the
18S protist and ITS fungal mock-community data.  Results from an additional BaNaNA 18S analysis are recorded in
`expected_genera_recovery.tsv`, `expected_genera_recovery_details.tsv`, and
`logs/banana/`, but are not plotted.

For a command-by-command guide beginning with acquisition of `ERR15092745`,
18S coordinate extraction, placement of the study-provided ITS reads, database
setup, the additional BaNaNA installation, partial targets, and audit of the
recorded run, see [`HOWTO.md`](HOWTO.md).

## Design

- Exact nested subsamples of 100, 500, 2,500, 12,500, 60,000, and 300,000 non-empty
  reads were generated with seed 727271. Consequently, every smaller dataset is
  a subset of the corresponding larger dataset.

- Savont 0.6.3 was run with its default ASV parameters. Its accepted read-length
  interval was set to 1,000–2,600 nt for 18S and 300–1,000 nt for ITS.

- UNOISE3 dereplication used `-minuniquesize 2`, followed by either the default
  minimum abundance or `-minsize 3`. Read-to-zOTU abundance mapping used a
  separate dereplication with `-minuniquesize 1`.

- Savont ASVs and UNOISE zOTUs were classified through Savont's
  minimap2-based classifier against SILVA 138.2 (18S) or UNITE 10.0 (ITS).

- A genus was counted as recovered when at least one reported sequence received
  the expected genus, a configured database alias, or a binomial species name
  beginning with that genus.

The analysis uses six expected genera: *Prymnesium*, *Euglena*,
*Chlorella*, *Paramecium*, *Gymnodinium*, and *Cryptomonas*.

## Results

At 300,000 reads, Savont recovered all 6/6 expected 18S genera and 9/10
expected ITS genera. The 18S sample is heavily dominated by *Prymnesium* and
contains several unexpected eukaryotic classifications, consistent with
contamination.

## Additional BaNaNA analysis

BaNaNA was also run on the 18S data using its Kit 14 workflow at pinned commit
`e5e65e6deb71eaeb20a9a4cd1161b3d2b64f02a2`. The input was already extracted
to the 18S region, so the redundant Barrnap and `extracting_rrna.py` stages were
bypassed. All downstream stages were retained. BaNaNA recovered 5/6 expected
genera at 300,000 reads. These results are retained for supplementary reporting
but are omitted from the plotted figure.

## Key outputs

- `expected_genera_recovery.tsv`: recovery count for every marker, method, and
  depth, including the additional BaNaNA results.
- `expected_genera_recovery_details.tsv`: per-expected-genus calls and matching
  database names.
- `figures/expected_genera_recovered.pdf` and `.svg`: supplementary figure
  (Savont vs. UNOISE3; BaNaNA is tabulated but not plotted).
- `logs/`: tool logs.

## Reproduction

From the repository root:

```bash
XDG_CACHE_HOME=/tmp/snakemake-cache \
MPLCONFIGDIR=/tmp/matplotlib-its-18s \
snakemake \
  -s reviewer_analyses/its_18s_benchmark/Snakefile \
  --cores 20 --rerun-incomplete
```

Executable names and input paths are configured in `config.yaml`. Place the
three external inputs listed in `input/input_manifest.tsv` under `input/`
before running. Raw reads, reference databases, the BaNaNA checkout and
environment, and full run trees are excluded from the repository.

The SILVA database can be obtained with:

```bash
savont download \
  --location reviewer_analyses/its_18s_benchmark/databases \
  --dbs silva-138.2
```

The workflow prepares the Savont-compatible UNITE 10.0 database from the exact
archive named in `input/input_manifest.tsv` after it is placed under `input/`.

The full workflow includes the additional BaNaNA analysis. `HOWTO.md` gives
commands to recreate `external/BaNaNA` at the pinned commit and `envs/banana`
from `banana_environment.yaml`.
