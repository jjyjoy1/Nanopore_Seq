rule extract_splice_junctions:
    input:
        bam = config["output"]["align_dir"] + "/{sample}.bam"
    output:
        junctions = config["output"]["splicing_dir"] + "/{sample}_junctions.csv"
    params:
        min_support = config["params"]["splicing"]["min_support"]
    log:
        config["output"]["splicing_dir"] + "/logs/{sample}_junctions.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/extract_junctions.py"

rule compare_splice_junctions:
    input:
        tumor_junctions = lambda wildcards: expand(
            config["output"]["splicing_dir"] + "/{sample}_junctions.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "tumor"]
        ),
        normal_junctions = lambda wildcards: expand(
            config["output"]["splicing_dir"] + "/{sample}_junctions.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "normal"]
        )
    output:
        tumor_specific = config["output"]["splicing_dir"] + "/tumor_specific_junctions.csv",
        differential = config["output"]["splicing_dir"] + "/differential_junctions.csv",
        plots_dir = directory(config["output"]["splicing_dir"] + "/plots")
    log:
        config["output"]["splicing_dir"] + "/logs/junction_comparison.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/splice_analysis.py"
