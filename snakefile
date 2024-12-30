# Snakemake pipeline for processing FASTQ files

import glob

# Automatically detect samples from FASTQ files
SAMPLES = [fastq.split(".fastq")[0] for fastq in glob.glob("*.fastq")]

rule all:
    input:
        expand("reports/fastqc/{sample}_fastqc.zip", sample=SAMPLES),
        expand("reports/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        expand("processed/{sample}_trimmed.fastq", sample=SAMPLES),
        expand("reports/kraken2/{sample}_kraken2_report.txt", sample=SAMPLES),
        expand("reports/amrfinder/{sample}_amrfinder_report.txt", sample=SAMPLES),
        expand("processed/{sample}_mapped.bam", sample=SAMPLES),
        expand("reports/qualimap/{sample}_qualimap_report", sample=SAMPLES),
        directory("resfinder_db/resfinder"),
        expand("reports/resfinder/{sample}_resfinder_report", sample=SAMPLES)

# Rule to download and build the Kraken2 database
rule kraken2_download_db:
    output:
        db_dir="kraken2_db",
        taxo="kraken2_db/taxo.k2d",
        hash="kraken2_db/hash.k2d",
        opts="kraken2_db/opts.k2d"
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        mkdir -p {output.db_dir} && \
        kraken2-build --standard --db {output.db_dir} && \
        kraken2-build --build --db {output.db_dir}
        """

# Rule to run FastQC
rule fastqc:
    input:
        fastq="{sample}.fastq"
    output:
        "reports/fastqc/{sample}_fastqc.zip",
        "reports/fastqc/{sample}_fastqc.html"
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc {input.fastq} -o reports/fastqc"

# Rule to run Porechop-ABI
rule porechop:
    input:
        fastq="{sample}.fastq"
    output:
        "processed/{sample}_trimmed.fastq"
    conda:
        "envs/porechop.yaml"
    shell:
        "porechop_abi --ab_initio -i {input.fastq} -o {output}"

# Rule to run Kraken2
rule kraken2:
    input:
        trimmed="processed/{sample}_trimmed.fastq",
        db="kraken2_db/taxo.k2d"
    output:
        "reports/kraken2/{sample}_kraken2_report.txt"
    conda:
        "envs/kraken2.yaml"
    shell:
        "kraken2 --db kraken2_db --output {output} {input.trimmed}"


# Rule to run ResFinder
rule resfinder:
    input:
        trimmed="processed/{sample}_trimmed.fastq",
        db="resfinder_db/resfinder/"
    output:
        directory("reports/resfinder/{sample}_resfinder_report")
    conda:
        "envs/resfinder.yaml"
    shell:
        """
        run_resfinder.py -ifq {input.trimmed} --nanopore -o reports/resfinder/{wildcards.sample}_resfinder_report --acquired -db_res {input.db}
        """

# Rule to run AMRFinder
rule amrfinder:
    input:
        trimmed="processed/{sample}_trimmed.fasta"
    output:
        "reports/amrfinder/{sample}_amrfinder_report.txt"
    conda:
        "envs/amrfinder.yaml"
    shell:
        """
        amrfinder -u
        amrfinder -n {input.trimmed} -o {output}
        """

# Rule to map reads using Minimap2
rule minimap2:
    input:
        trimmed="processed/{sample}_trimmed.fastq",
        reference="reference.fasta"
    output:
        bam="processed/{sample}_mapped.bam"
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        minimap2 -a {input.reference} {input.trimmed} | samtools view -Sb - > {output.bam}.tmp
        samtools sort -o {output.bam} -m 20G --output-fmt BAM --reference {input.reference} {output.bam}.tmp
        rm -v {output.bam}.tmp
        """

# Rule to run Qualimap
rule qualimap:
    input:
        bam="processed/{sample}_mapped.bam"
    output:
        directory("reports/qualimap/{sample}_qualimap_report")
    conda:
        "envs/qualimap.yaml"
    shell:
        "qualimap bamqc -bam {input.bam} -outdir reports/qualimap/{wildcards.sample}_qualimap_report"

# Rule to convert fastq to fasta
rule fastq_to_fasta:
    input:
        fastq="processed/{sample}_trimmed.fastq"
    output:
        fasta="processed/{sample}_trimmed.fasta"
    conda:
        "envs/seqtk.yaml"
    shell:
        """
        seqtk seq -A {input.fastq} > {output.fasta}
        """

# Rule to download ResFinder database
rule resfinder_download_db:
    output:
        directory("resfinder_db/resfinder")
    conda:
        "envs/resfinder.yaml"
    shell:
        """
        mkdir -p {output} && git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git {output} && cd {output}
        python3 INSTALL.py `which kma` non_interactive
        """
