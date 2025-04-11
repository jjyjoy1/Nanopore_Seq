rule qc_raw_reads:
    input:
        fastq = lambda wildcards: config["samples"][wildcards.sample]["fastq"]
    output:
        html = config["output"]["qc_dir"] + "/{sample}_nanoplot.html",
        stats = config["output"]["qc_dir"] + "/{sample}_nanoplot_stats.txt"
    log:
        config["output"]["qc_dir"] + "/logs/{sample}_nanoplot.log"
    conda:
        "../envs/base.yaml"
    shell:
        """
        NanoPlot --fastq {input.fastq} \
            --outdir $(dirname {output.html}) \
            --prefix {wildcards.sample}_nanoplot \
            2> {log}
        """

rule qc_aligned_reads:
    input:
        bam = config["output"]["align_dir"] + "/{sample}.bam"
    output:
        metrics = config["output"]["qc_dir"] + "/{sample}_alignment_metrics.txt"
    log:
        config["output"]["qc_dir"] + "/logs/{sample}_alignment_qc.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/qc_metrics.py"

rule multiqc:
    input:
        expand(config["output"]["qc_dir"] + "/{sample}_nanoplot_stats.txt", sample=SAMPLES),
        expand(config["output"]["qc_dir"] + "/{sample}_alignment_metrics.txt", sample=SAMPLES)
    output:
        report = config["output"]["qc_dir"] + "/multiqc_report.html"
    log:
        config["output"]["qc_dir"] + "/logs/multiqc.log"
    conda:
        "../envs/base.yaml"
    shell:
        """
        multiqc {config[output][qc_dir]} \
            -o {config[output][qc_dir]} \
            -f \
            -n multiqc_report.html \
            2> {log}
        """

