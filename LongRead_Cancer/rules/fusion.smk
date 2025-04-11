rule extract_gene_positions:
    input:
        annotation = config["reference"]["annotation"]
    output:
        gene_positions = config["output"]["fusion_dir"] + "/gene_positions.pkl"
    log:
        config["output"]["fusion_dir"] + "/logs/gene_positions.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/extract_gene_positions.py"

rule detect_fusion_candidates:
    input:
        bam = config["output"]["align_dir"] + "/{sample}.bam",
        gene_positions = config["output"]["fusion_dir"] + "/gene_positions.pkl"
    output:
        fusions = config["output"]["fusion_dir"] + "/{sample}_fusions.csv"
    params:
        min_support = config["params"]["fusion"]["min_support"]
    log:
        config["output"]["fusion_dir"] + "/logs/{sample}_fusion_detection.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/fusion_detection.py"

rule filter_tumor_specific_fusions:
    input:
        tumor_fusions = lambda wildcards: expand(
            config["output"]["fusion_dir"] + "/{sample}_fusions.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "tumor"]
        ),
        normal_fusions = lambda wildcards: expand(
            config["output"]["fusion_dir"] + "/{sample}_fusions.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "normal"]
        )
    output:
        tumor_specific = config["output"]["fusion_dir"] + "/tumor_specific_fusions.csv",
        plots_dir = directory(config["output"]["fusion_dir"] + "/plots")
    log:
        config["output"]["fusion_dir"] + "/logs/fusion_filtering.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/filter_fusions.py"


