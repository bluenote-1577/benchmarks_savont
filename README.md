# Savont Benchmarking

Benchmarking pipelines for **savont**, an ASV generation method for long-read 16S rRNA sequencing (ONT and PacBio HiFi).

## Data Acquisition

### Read Files

Download from SRA and place in the appropriate directories:

**Mock Community (Zymo):**

| File Path | SRA Accession |
|-----------|---------------|
| `reads/morten-zymo-simple/v1v8_zmock1.fastq.gz` | [SRR36567714](https://www.ncbi.nlm.nih.gov/sra/SRR36567714) |
| `reads/morten-zymo-simple-opr/v1v8_zmock1_opr.fastq.gz` | [SRR36567737](https://www.ncbi.nlm.nih.gov/sra/SRR36567737) |
| `reads/ziels_pb_zymo/ziels_pb_zymo.fastq.gz` | [ERR3813246](https://www.ebi.ac.uk/ena/browser/view/ERR3813246) |
| `reads/ziels_pb_zymo/ziels_pb_zymo_opr.fastq.gz` | [ERR3813246](https://www.ebi.ac.uk/ena/browser/view/ERR3813246) (trimmed) |

**Real Environmental Samples (ONT):**

| File Path | SRA Accession | Sample Type |
|-----------|---------------|-------------|
| `reads/real_data/AD1_ont.fastq.gz` | [SRR36567675](https://www.ncbi.nlm.nih.gov/sra/SRR36567675) | Anaerobic Digester |
| `reads/real_data/WWTP1_ont.fastq.gz` | [SRR36567683](https://www.ncbi.nlm.nih.gov/sra/SRR36567683) | Wastewater |
| `reads/real_data/Zfecal1_ont.fastq.gz` | [SRR36567672](https://www.ncbi.nlm.nih.gov/sra/SRR36567672) | Fecal |
| `reads/real_data/Soil1_ont.fastq.gz` | [SRR36567678](https://www.ncbi.nlm.nih.gov/sra/SRR36567678) | Soil |

**Note**: For HiFi V1V8 data, this is obtained by trimming the V8 region onward using cutadapt, so it's the same file as the operon reads.

### Reference Data

Included in `morten_data/`:
- `V1V8_bacterial_amplicons.fasta.gz` — V1-V8 reference amplicons (D6300 Zymo)
- `OPR_bacterial_amplicons.fasta.gz` — Operon reference amplicons

## Environment Setup

```bash
# Base environment (savont, minimap2, samtools, seqtk, porechop, etc.)
conda env create -f environments/base.yaml -n savont_bench

# Cutadapt environment
conda env create -f environments/cutadapt.yaml -n cutadapt

# QIIME2 environment (for DADA2 on HiFi data)
conda env create -f environments/qiime2-amplicon-2024.10.yaml -n qiime2-amplicon-2024.10
```

### Required Software

- **savont** — install from the savont repository
- **usearch** — required for UNOISE3 (free 32-bit or licensed 64-bit)
- **porechop**, **cutadapt**, **seqtk**, **minimap2**, **samtools** — via conda

## Running Benchmarks

| Config | Description |
|--------|-------------|
| `config_fl.yaml` | Full-length V1-V8 amplicons |
| `config_opr.yaml` | rRNA operon reads |
| `config_real.yaml` | Real environmental samples |

```bash
conda activate savont_bench

# Mock community benchmarks
snakemake -j 20 --use-conda ##change config files manually in Snakefile (fl or opr)

# Real data benchmark
snakemake -s snakefile_real.snk -j 20 --use-conda
```

The pipeline runs: adapter trimming (porechop) → primer trimming (cutadapt) → subsampling (seqtk) → ASV generation (savont, UNOISE3, DADA2) → alignment to reference (minimap2).

Results are written to `results/` with timing benchmarks in `results/benchmarks/`.

## Analysis Notebooks

| Notebook | Description |
|----------|-------------|
| `analyze_paf.ipynb` | ASV detection accuracy, false positives, runtime/memory benchmarks |
| `real_notebook/real_analysis.ipynb` | Real data ASV counts and mapping fidelity |
| `emu_analysis_real/divergence_analysis.ipynb` | EMU divergence vs abundance analysis |

**Notes:**
- `analyze_paf.ipynb` can be run directly with pre-computed results in `results/`
- `real_notebook/real_analysis.ipynb` requires running `snakemake -s snakefile_real.snk` first
- For the EMU notebook, first run `emu abundance` with `--keep-files --keep-counts`, then run `emu_analysis.py`

## Reviewer analyses

Compact final tables and plotting code for the revision experiments are in
`reviewer_analyses/`. They cover abundance-stratified HiFi agreement, minimum
cluster size, primary clustering identity, simulated read accuracy, HiFi zOTU
recall by observed ONT depth, and ITS/18S expected-genus recovery. Large raw and
temporary outputs are intentionally omitted.

To regenerate all six reviewer figures from the included tables:

```bash
conda env create -f reviewer_analyses/plotting_environment.yaml
conda activate savont-reviewer-plots
MPLCONFIGDIR=/tmp/matplotlib-reviewer \
python reviewer_analyses/plot_reviewer_figures.py
```

See `reviewer_analyses/README.md` for the input table associated with each
figure and commands for regenerating selected panels. A command-by-command
guide for rebuilding the ITS/18S inputs and full benchmark is available in
`reviewer_analyses/its_18s_benchmark/HOWTO.md`; compact logs from the actual
reviewer-analysis runs are retained alongside their result tables.

## Directory Structure

```
├── Snakefile                    # Main pipeline (mock community)
├── snakefile_real.snk           # Pipeline for real samples
├── config_*.yaml                # Configuration files
├── environments/                # Conda environment files
├── reads/                       # Input reads (download from SRA)
├── morten_data/                 # Reference sequences
├── morten-zotu-data/            # Pre-computed ZOTU data
├── results/                     # Pipeline outputs
├── analyze_paf.ipynb            # Main analysis notebook
├── real_notebook/               # Real data analysis
├── reviewer_analyses/           # Compact revision tables and plotting code
└── emu_analysis_real/           # EMU analysis
```

### AI usage

We used Claude Sonnet v4.6 and GPT-Sol 5.6 to help with scripting and organize this repository. 
