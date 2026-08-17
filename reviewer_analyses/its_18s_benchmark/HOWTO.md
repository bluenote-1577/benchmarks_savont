# HOWTO: reproduce the ITS and 18S benchmark

## 1. Install the general tools

The workflow needs Python with pandas/matplotlib/seaborn, Snakemake, Savont,
USEARCH, minimap2, HMMER/Barrnap, and seqtk. USEARCH must be obtained under its
own license and made available as `usearch`. Install Savont separately and make
`savont` available on `PATH`.

The main repository environment contains most supporting programs:

```bash
conda env create -f environments/base.yaml -n savont_bench
conda activate savont_bench
```

The recorded manuscript versions are in `banana_provenance.tsv`, the
analysis-specific `software_version*.txt` files, and the repository's software
version table.

## 2. Reconstruct the extracted 18S FASTQ

The public starting run is ENA/SRA accession `ERR15092745` from project
`PRJEB89945`. Download this accession with an SRA/ENA downloader.

The analyzed file was not the entire 18S–D2/28S amplicon. It was the 18S
region extracted before benchmarking. The extraction used Barrnap 0.9's
eukaryotic `18S_rRNA` model (`RF01960`) directly with Barrnap's nhmmer engine,
then sliced the reported coordinates from FASTQ while preserving qualities and
reverse-complementing reverse-strand hits. Direct nhmmer invocation was used so
the exact hit coordinates were available.

```bash
seqtk seq -A reviewer_analyses/its_18s_benchmark/input/ERR15092745.fastq \
  > reviewer_analyses/its_18s_benchmark/input/ERR15092745.fasta

BARRNAP_EXE=$(readlink -f "$(command -v barrnap)")
BARRNAP_EUK_HMM="$(dirname "$BARRNAP_EXE")/../db/euk.hmm"
hmmfetch "$BARRNAP_EUK_HMM" 18S_rRNA \
  > reviewer_analyses/its_18s_benchmark/input/18S_rRNA.hmm

mkdir -p reviewer_analyses/its_18s_benchmark/logs/input_preparation
nhmmer --noali --cpu 32 -E 1e-10 \
  --tblout reviewer_analyses/its_18s_benchmark/input/ERR15092745.18S.tbl \
  reviewer_analyses/its_18s_benchmark/input/18S_rRNA.hmm \
  reviewer_analyses/its_18s_benchmark/input/ERR15092745.fasta \
  > reviewer_analyses/its_18s_benchmark/logs/input_preparation/ERR15092745.18S.nhmmer.log

python reviewer_analyses/its_18s_benchmark/scripts/extract_18s_hits.py \
  --tblout reviewer_analyses/its_18s_benchmark/input/ERR15092745.18S.tbl \
  --fastq reviewer_analyses/its_18s_benchmark/input/ERR15092745.fastq \
  --output reviewer_analyses/its_18s_benchmark/input/ERR15092745.18S.fastq
```

The expected result contains 580,918 reads.

This extraction is why the BaNaNA comparison bypasses only its Barrnap and
`extracting_rrna.py` stages. Running those stages again would extract an
already-extracted sequence a second time. All downstream BaNaNA stages remain
enabled in `scripts/run_banana_otu.py`.

## 3. Place the trimmed ITS FASTQ

The ITS input is the study-provided, already primer-trimmed and cut ONT file,
not a file derived by this repository. Copy it to:

```text
reviewer_analyses/its_18s_benchmark/input/T1_B2_all_trimmed.fq
```

The exact file has 855,782 FASTQ records, including three zero-length records.

The upstream primer-trimming inputs and primer sequences were not available in
this repository, so the workflow begins from this trimmed FASTQ.

## 4. Prepare SILVA, UNITE, BaNaNA, and PR2

Download SILVA 138.2 in Savont format:

```bash
savont download \
  --location reviewer_analyses/its_18s_benchmark/databases \
  --dbs silva-138.2
```

The Snakemake rule `prepare_unites` converts its
long pipe-delimited headers into Savont's FASTA plus taxonomy-table format.

For BaNaNA, clone the exact source commit and create the pinned environment:

```bash
git clone https://github.com/ibe-uw/BaNaNA.git \
  reviewer_analyses/its_18s_benchmark/external/BaNaNA
git -C reviewer_analyses/its_18s_benchmark/external/BaNaNA \
  checkout e5e65e6deb71eaeb20a9a4cd1161b3d2b64f02a2

conda env create \
  -p reviewer_analyses/its_18s_benchmark/envs/banana \
  -f reviewer_analyses/its_18s_benchmark/banana_environment.yaml
```

Download PR2 5.1.1's SSU UTAX FASTA from the PR2 release archive, decompress
it, and place it at:

```text
reviewer_analyses/its_18s_benchmark/databases/pr2-5.1.1/pr2_version_5.1.1_SSU_UTAX.fasta
```

PR2 is used only for BaNaNA's reference-chimera filtering. All final 18S
taxonomic scoring uses SILVA so methods are compared through one classifier.


## 5. Run the workflow

Check the DAG first:

```bash
snakemake \
  -s reviewer_analyses/its_18s_benchmark/Snakefile \
  --cores 1 --dry-run
```

Run all depths and methods:

```bash
XDG_CACHE_HOME=/tmp/snakemake-cache \
MPLCONFIGDIR=/tmp/matplotlib-its-18s \
snakemake \
  -s reviewer_analyses/its_18s_benchmark/Snakefile \
  --cores 20 --rerun-incomplete --printshellcmds
```


## 6. Audit the manuscript run

- `logs/savont/`: all 12 Savont ASV logs.
- `logs/unoise/`: dereplication, denoising, and abundance-mapping logs.
  A zero-byte OTU-table log is expected when its paired denoising log records
  zero zOTUs; in that case no USEARCH mapping command was invoked.
- `logs/banana/*.log.gz`: complete BaNaNA command/stage logs, compressed
  losslessly because the 300,000-read log is verbose. Inspect with `zless`.
- `logs/classification/`: harmonized Savont classification logs for every
  marker/method/depth combination.
- `logs/database/`: minimap2 index construction logs.
- `results/subsampled/*.audit.tsv`: source counts and length distributions.
- `results/banana/18s_*/run_summary.tsv`: read/cluster/consensus/chimera/OTU
  counts at each BaNaNA depth.
- `expected_genera_recovery_details.tsv`: the exact species/genus strings that
  caused each expected-genus recovery call.

These files demonstrate execution and make the scoring decisions auditable
without committing raw reads, reference databases, or multi-gigabyte temporary
run trees.
