rule extract_variants:
    input:
        bam = config["output"]["align_dir"] + "/{sample}.bam",
        vcf = config["reference"]["variants"]
    output:
        variants = config["output"]["ase_dir"] + "/{sample}_read_variants.pkl"
    log:
        config["output"]["ase_dir"] + "/logs/{sample}_extract_variants.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/extract_variants.py"

rule phase_variants:
    input:
        variants = config["output"]["ase_dir"] + "/{sample}_read_variants.pkl"
    output:
        phased = config["output"]["ase_dir"] + "/{sample}_phased_blocks.pkl"
    log:
        config["output"]["ase_dir"] + "/logs/{sample}_phase_variants.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/phase_variants.py"

rule quantify_ase:
    input:
        bam = config["output"]["align_dir"] + "/{sample}.bam",
        phased = config["output"]["ase_dir"] + "/{sample}_phased_blocks.pkl",
        annotation = config["reference"]["annotation"]
    output:
        ase_results = config["output"]["ase_dir"] + "/{sample}_ase_results.csv",
        plots_dir = directory(config["output"]["ase_dir"] + "/{sample}_ase_plots")
    params:
        min_coverage = config["params"]["ase"]["min_coverage"]
    log:
        config["output"]["ase_dir"] + "/logs/{sample}_quantify_ase.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/quantify_ase.py"

