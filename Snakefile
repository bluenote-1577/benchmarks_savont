configfile: "config_opr.yaml"

import math

# Wildcard constraints
wildcard_constraints:
    sample="[^/]+",
    fraction="[0-9.]+",
    method="(savont|unoise|unoise_min3)"

# Get all samples, fractions, and methods
SAMPLE_META = config["samples"]
SAMPLES     = list(SAMPLE_META.keys())
ONT_SAMPLES = [s for s, m in SAMPLE_META.items() if m["platform"] == "ont"]
PB_SAMPLES  = [s for s, m in SAMPLE_META.items() if m["platform"] == "pb"]
METHODS     = ["savont", "unoise", "unoise_min3"]
#METHODS = ["unoise", "unoise_min3"]

# All valid (sample, fraction) pairs — fractions differ per sample
SAMPLE_FRACTIONS = [
    (s, f)
    for s, m in SAMPLE_META.items()
    for f in m["fractions"]
]

# Final outputs
rule all:
    input:
        [f"results/alignments/{s}_frac{f}_{method}.paf"
         for s, f in SAMPLE_FRACTIONS
         for method in METHODS],
        [f"results/alignments/{s}_frac{f}_dada2.paf"
         for s, f in SAMPLE_FRACTIONS
         if s in PB_SAMPLES],
        [f"results/coverage/{s}_frac{f}.tsv"
         for s, f in SAMPLE_FRACTIONS]

# Step 1: Run porechop on raw reads
rule porechop:
    input:
        lambda wildcards: SAMPLE_META[wildcards.sample]["reads"]
    output:
        temp("results/porechop/{sample}.fastq")
    threads:
        config["threads"]
    log:
        "logs/porechop/{sample}.log"
    conda:
        "savont_bench"
    shell:
        """
        porechop -i {input} -o {output} -t {threads} 2>&1 | tee {log}
        """

# Step 2: Run cutadapt to trim primers
rule cutadapt:
    input:
        "results/porechop/{sample}.fastq"
    output:
        "results/cutadapt/{sample}.fastq"
    params:
        primer=config["primers"][config["primer_type"]],
        error_rate=config["cutadapt_error_rate"]
    threads:
        config["threads"]
    log:
        "logs/cutadapt/{sample}.log"
    conda:
        "cutadapt"
    shell:
        """
        cutadapt -g '{params.primer}' \
            --rc \
            {input} \
            -j {threads} \
            --discard-untrimmed \
            -e {params.error_rate} \
            -o {output} \
            2>&1 | tee {log}
        """

# Step 3: Subsample reads at different fractions
rule subsample:
    input:
        "results/cutadapt/{sample}.fastq"
    output:
        "results/subsampled/{sample}_frac{fraction}.fastq"
    params:
        seed=42
    log:
        "logs/subsample/{sample}_frac{fraction}.log"
    run:
        fraction = float(wildcards.fraction)
        if fraction == 1.0:
            # No subsampling, just copy
            shell("cp {input} {output}")
        else:
            # Use seqtk for subsampling
            shell(f"seqtk sample -s {{params.seed}} {{input}} {fraction} > {{output}} 2> {{log}}")

# Step 4a: Run UNOISE3 pipeline
rule unoise_uniques:
    input:
        "results/subsampled/{sample}_frac{fraction}.fastq"
    output:
        "results/unoise/{sample}_frac{fraction}/uniques.fa"
    params:
        min_size=config["min_unique_size"]
    log:
        "logs/unoise/{sample}_frac{fraction}_uniques.log"
    benchmark:
        "results/benchmarks/unoise_uniques/{sample}_frac{fraction}.tsv"
    shell:
        """
        usearch -fastx_uniques {input} \
            -sizeout \
            -relabel Uniq \
            -minuniquesize {params.min_size} \
            -fastaout {output} \
            2>&1 | tee {log} || true
        """

rule unoise_denoise:
    input:
        "results/unoise/{sample}_frac{fraction}/uniques.fa"
    output:
        "results/unoise/{sample}_frac{fraction}/zotus.fa"
    benchmark:
        "results/benchmarks/unoise_denoise/{sample}_frac{fraction}.tsv"
    log:
        "logs/unoise/{sample}_frac{fraction}_denoise.log"
    conda:
        "savont_bench"
    shell:
        """
        touch {output}
        usearch -unoise3 {input} \
            -zotus {output} \
            2>&1 | tee {log} || true
        """

# Step 4a-alt: Run UNOISE3 with minsize 3 (reuses uniques from standard unoise)
rule unoise_denoise_min3:
    input:
        "results/unoise/{sample}_frac{fraction}/uniques.fa"
    output:
        "results/unoise_min3/{sample}_frac{fraction}/zotus.fa"
    benchmark:
        "results/benchmarks/unoise_denoise_min3/{sample}_frac{fraction}.tsv"
    log:
        "logs/unoise_min3/{sample}_frac{fraction}_denoise.log"
    conda:
        "savont_bench"
    shell:
        """
        touch {output}
        usearch -unoise3 {input} \
            -zotus {output} \
            -minsize 3 \
            2>&1 | tee {log} || true
        """

# Step 4b: Run savont
rule savont:
    input:
        "results/subsampled/{sample}_frac{fraction}.fastq"
    output:
        "results/savont/{sample}_frac{fraction}/final_asvs.fasta"
    benchmark:
        "results/benchmarks/savont/{sample}_frac{fraction}.tsv"
    params:
        outdir="results/savont/{sample}_frac{fraction}",
        opr=config['primer_type'],
        hifi_flag=lambda wildcards: "--hifi" if wildcards.sample in PB_SAMPLES else ""
    threads:
        config["threads"]
    log:
        "logs/savont/{sample}_frac{fraction}.log"
    conda:
        "savont_bench"
    shell:
        """
        if [ "opr" = "{params.opr}" ]; then
          savont asv -t {threads} {input} --rrna-operon {params.hifi_flag} -o {params.outdir} 2>&1 | tee {log}
        else
          savont asv -t {threads} {input} {params.hifi_flag} -o {params.outdir} 2>&1 | tee {log}
        fi
        """

# Step 5: Run minimap2 alignment for unoise
rule minimap2_unoise:
    input:
        ref=config["reference"],
        query="results/unoise/{sample}_frac{fraction}/zotus.fa"
    output:
        "results/alignments/{sample}_frac{fraction}_unoise.paf"
    threads:
        config["threads"]
    log:
        "logs/minimap2/{sample}_frac{fraction}_unoise.log"
    conda:
        "savont_bench"
    shell:
        """
        minimap2 --cs --secondary=no \
            -t {threads} \
            {input.ref} \
            {input.query} \
            > {output} \
            2> {log}
        """

# Step 5: Run minimap2 alignment for savont
rule minimap2_savont:
    input:
        ref=config["reference"],
        query="results/savont/{sample}_frac{fraction}/final_asvs.fasta"
    output:
        "results/alignments/{sample}_frac{fraction}_savont.paf"
    threads:
        config["threads"]
    log:
        "logs/minimap2/{sample}_frac{fraction}_savont.log"
    conda:
        "savont_bench"
    shell:
        """
        minimap2 --cs --secondary=no \
            -t {threads} \
            {input.ref} \
            {input.query} \
            > {output} \
            2> {log}
        """

# Step 5: Run minimap2 alignment for unoise_min3
rule minimap2_unoise_min3:
    input:
        ref=config["reference"],
        query="results/unoise_min3/{sample}_frac{fraction}/zotus.fa"
    output:
        "results/alignments/{sample}_frac{fraction}_unoise_min3.paf"
    threads:
        config["threads"]
    log:
        "logs/minimap2/{sample}_frac{fraction}_unoise_min3.log"
    conda:
        "savont_bench"
    shell:
        """
        minimap2 --cs --secondary=no \
            -t {threads} \
            {input.ref} \
            {input.query} \
            > {output} \
            2> {log}
        """


# ── Read coverage against reference ─────────────────────────────────────────

rule read_coverage:
    input:
        ref=config["reference"],
        reads="results/subsampled/{sample}_frac{fraction}.fastq"
    output:
        "results/coverage/{sample}_frac{fraction}.tsv"
    threads:
        config["threads"]
    log:
        "logs/coverage/{sample}_frac{fraction}.log"
    conda:
        "savont_bench"
    shell:
        """
        minimap2 -ax lr:hq \
            -t {threads} \
            {input.ref} \
            {input.reads} \
            2>> {log} \
        | samtools sort -@ {threads} -o /tmp/{wildcards.sample}_frac{wildcards.fraction}.bam
        samtools index /tmp/{wildcards.sample}_frac{wildcards.fraction}.bam
        samtools coverage /tmp/{wildcards.sample}_frac{wildcards.fraction}.bam \
            2>> {log} \
        | sed '1s/^#//' \
        > {output}
        rm /tmp/{wildcards.sample}_frac{wildcards.fraction}.bam \
           /tmp/{wildcards.sample}_frac{wildcards.fraction}.bam.bai
        """


# ── DADA2 pipeline (PacBio CCS reads only) ──────────────────────────────────

# Step D1: Create QIIME2 manifest TSV pointing to subsampled reads
rule dada2_manifest:
    input:
        "results/subsampled/{sample}_frac{fraction}.fastq"
    output:
        "results/dada2/{sample}_frac{fraction}/manifest.tsv"
    run:
        import os
        abs_path = os.path.abspath(str(input[0]))
        with open(str(output[0]), "w") as fh:
            fh.write("sample-id\tabsolute-filepath\n")
            fh.write(f"{wildcards.sample}\t{abs_path}\n")

# Step D2: Import reads into QIIME2 artifact
rule dada2_import:
    input:
        "results/dada2/{sample}_frac{fraction}/manifest.tsv"
    output:
        "results/dada2/{sample}_frac{fraction}/sequences.qza"
    log:
        "logs/dada2/{sample}_frac{fraction}_import.log"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        """
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path {input} \
            --input-format SingleEndFastqManifestPhred33V2 \
            --output-path {output} \
            2>&1 | tee {log}
        """

# Step D3: Denoise with DADA2 CCS
rule dada2_denoise:
    input:
        "results/dada2/{sample}_frac{fraction}/sequences.qza"
    output:
        rep_seqs="results/dada2/{sample}_frac{fraction}/denoise/representative_sequences.qza",
        table="results/dada2/{sample}_frac{fraction}/denoise/table.qza",
        stats="results/dada2/{sample}_frac{fraction}/denoise/denoising_stats.qza"
    benchmark:
        "results/benchmarks/dada2_denoise/{sample}_frac{fraction}.tsv"
    params:
        outdir="results/dada2/{sample}_frac{fraction}/denoise"
    threads:
        config["threads"]
    log:
        "logs/dada2/{sample}_frac{fraction}_denoise.log"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
      ### NOTE: using a modified version of qiime dada2 denoise-ccs that allows no primers to be specified. 
        """
        rm -r {params.outdir} || true
        qiime dada2 denoise-ccs \
            --i-demultiplexed-seqs {input} \
            --output-dir {params.outdir} \
            --p-n-threads {threads} \
            2>&1 | tee {log}
        """

# Step D4: Export representative sequences to FASTA
rule dada2_export:
    input:
        "results/dada2/{sample}_frac{fraction}/denoise/representative_sequences.qza"
    output:
        "results/dada2/{sample}_frac{fraction}/exported/dna-sequences.fasta"
    params:
        outdir="results/dada2/{sample}_frac{fraction}/exported"
    log:
        "logs/dada2/{sample}_frac{fraction}_export.log"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {params.outdir} \
            2>&1 | tee {log}
        """

# Step D5: Align DADA2 ASVs to reference
rule minimap2_dada2:
    input:
        ref=config["reference"],
        query="results/dada2/{sample}_frac{fraction}/exported/dna-sequences.fasta"
    output:
        "results/alignments/{sample}_frac{fraction}_dada2.paf"
    threads:
        config["threads"]
    log:
        "logs/minimap2/{sample}_frac{fraction}_dada2.log"
    conda:
        "savont_bench"
    shell:
        """
        minimap2 --cs --secondary=no \
            -t {threads} \
            {input.ref} \
            {input.query} \
            > {output} \
            2> {log}
        """
