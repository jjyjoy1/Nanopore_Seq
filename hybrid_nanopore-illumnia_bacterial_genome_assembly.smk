# Snakefile for hybrid bacterial genome assembly using Nanopore and Illumina data

# Configuration file (save as config.yaml)
configfile: "config.yaml"

# Define final output files
rule all:
    input:
        final_assembly = "results/assembly/{sample}/final_assembly.fasta",
        quast_report = "results/qc/{sample}/quast/report.html",
        busco_report = "results/qc/{sample}/busco/short_summary.txt"

# Quality control for Nanopore reads
rule nanopore_qc:
    input:
        nanopore_reads = lambda wildcards: config["samples"][wildcards.sample]["nanopore"]
    output:
        qc_report = "results/qc/{sample}/nanoplot/NanoPlot-report.html",
        filtered_reads = "results/cleaned_reads/{sample}/nanopore_filtered.fastq.gz"
    log:
        "logs/nanopore_qc/{sample}.log"
    threads: 8
    shell:
        """
        # Generate QC report with NanoPlot
        mkdir -p results/qc/{wildcards.sample}/nanoplot
        NanoPlot --fastq {input.nanopore_reads} -o results/qc/{wildcards.sample}/nanoplot -t {threads} 2> {log}
        
        # Filter reads with Filtlong
        mkdir -p results/cleaned_reads/{wildcards.sample}
        filtlong --min_length 1000 --keep_percent 90 {input.nanopore_reads} | gzip > {output.filtered_reads} 2>> {log}
        """

# Quality control and adapter trimming for Illumina reads
rule illumina_qc:
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["illumina_r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["illumina_r2"]
    output:
        r1_trimmed = "results/cleaned_reads/{sample}/illumina_R1_trimmed.fastq.gz",
        r2_trimmed = "results/cleaned_reads/{sample}/illumina_R2_trimmed.fastq.gz",
        html = "results/qc/{sample}/fastp/fastp.html",
        json = "results/qc/{sample}/fastp/fastp.json"
    log:
        "logs/illumina_qc/{sample}.log"
    threads: 8
    shell:
        """
        mkdir -p results/qc/{wildcards.sample}/fastp
        fastp -i {input.r1} -I {input.r2} \
            -o {output.r1_trimmed} -O {output.r2_trimmed} \
            --html {output.html} --json {output.json} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --unqualified_percent_limit 40 \
            --n_base_limit 5 \
            --length_required 50 \
            --thread {threads} \
            2> {log}
        """

# De novo assembly with Flye using Nanopore reads
rule flye_assembly:
    input:
        nanopore_reads = "results/cleaned_reads/{sample}/nanopore_filtered.fastq.gz"
    output:
        assembly = "results/assembly/{sample}/flye/assembly.fasta",
        info = "results/assembly/{sample}/flye/assembly_info.txt"
    log:
        "logs/flye/{sample}.log"
    params:
        genome_size = lambda wildcards: config["samples"][wildcards.sample]["genome_size"],
        out_dir = "results/assembly/{sample}/flye"
    threads: 16
    shell:
        """
        flye --nano-raw {input.nanopore_reads} \
            --genome-size {params.genome_size} \
            --out-dir {params.out_dir} \
            --threads {threads} \
            2> {log}
        """

# Polish the assembly with Illumina reads using Pilon
rule polish_with_illumina:
    input:
        assembly = "results/assembly/{sample}/flye/assembly.fasta",
        r1_trimmed = "results/cleaned_reads/{sample}/illumina_R1_trimmed.fastq.gz",
        r2_trimmed = "results/cleaned_reads/{sample}/illumina_R2_trimmed.fastq.gz"
    output:
        polished_assembly = "results/assembly/{sample}/pilon/pilon.fasta",
        bam = temp("results/assembly/{sample}/pilon/aligned.bam"),
        bai = temp("results/assembly/{sample}/pilon/aligned.bai")
    log:
        "logs/pilon/{sample}.log"
    threads: 16
    shell:
        """
        # Index the assembly
        mkdir -p results/assembly/{wildcards.sample}/pilon
        bwa index {input.assembly}
        
        # Align Illumina reads to the assembly
        bwa mem -t {threads} {input.assembly} {input.r1_trimmed} {input.r2_trimmed} | \
            samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        
        # Polish with Pilon
        pilon --genome {input.assembly} \
            --frags {output.bam} \
            --output pilon \
            --outdir results/assembly/{wildcards.sample}/pilon \
            --changes \
            --fix bases \
            --threads {threads} \
            2> {log}
        """

# Further polishing with Racon and medaka
rule polish_with_nanopore:
    input:
        assembly = "results/assembly/{sample}/pilon/pilon.fasta",
        nanopore_reads = "results/cleaned_reads/{sample}/nanopore_filtered.fastq.gz"
    output:
        racon_assembly = temp("results/assembly/{sample}/racon/racon.fasta"),
        final_assembly = "results/assembly/{sample}/final_assembly.fasta"
    log:
        "logs/polish_nanopore/{sample}.log"
    threads: 16
    shell:
        """
        # Align Nanopore reads to the Pilon-polished assembly
        mkdir -p results/assembly/{wildcards.sample}/racon
        minimap2 -ax map-ont -t {threads} {input.assembly} {input.nanopore_reads} > results/assembly/{wildcards.sample}/racon/aligned.sam 2> {log}
        
        # Polish with Racon
        racon -t {threads} {input.nanopore_reads} results/assembly/{wildcards.sample}/racon/aligned.sam {input.assembly} > {output.racon_assembly} 2>> {log}
        
        # Polish with Medaka
        mkdir -p results/assembly/{wildcards.sample}/medaka
        medaka_consensus -i {input.nanopore_reads} -d {output.racon_assembly} -o results/assembly/{wildcards.sample}/medaka -t {threads} 2>> {log}
        
        # Copy the final assembly
        cp results/assembly/{wildcards.sample}/medaka/consensus.fasta {output.final_assembly}
        """

# Quality assessment with QUAST
rule quast:
    input:
        assembly = "results/assembly/{sample}/final_assembly.fasta",
        r1 = "results/cleaned_reads/{sample}/illumina_R1_trimmed.fastq.gz",
        r2 = "results/cleaned_reads/{sample}/illumina_R2_trimmed.fastq.gz"
    output:
        report = "results/qc/{sample}/quast/report.html"
    log:
        "logs/quast/{sample}.log"
    params:
        outdir = "results/qc/{sample}/quast"
    threads: 8
    shell:
        """
        quast.py {input.assembly} \
            -o {params.outdir} \
            --threads {threads} \
            --pe1 {input.r1} \
            --pe2 {input.r2} \
            2> {log}
        """

# Assess completeness with BUSCO
rule busco:
    input:
        assembly = "results/assembly/{sample}/final_assembly.fasta"
    output:
        summary = "results/qc/{sample}/busco/short_summary.txt"
    log:
        "logs/busco/{sample}.log"
    params:
        outdir = "results/qc/{sample}/busco",
        lineage = lambda wildcards: config["samples"][wildcards.sample]["busco_lineage"]
    threads: 8
    shell:
        """
        busco -i {input.assembly} \
            -o busco \
            -l {params.lineage} \
            -m genome \
            -c {threads} \
            --out_path {params.outdir} \
            2> {log}
        
        # Copy the summary file to the expected location
        cp {params.outdir}/busco/short_summary*.txt {output.summary}
        """
