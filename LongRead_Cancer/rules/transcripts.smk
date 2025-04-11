rule identify_transcripts:
    input:
        bam = config["output"]["align_dir"] + "/{sample}.bam",
        reference = config["reference"]["genome"],
        annotation = config["reference"]["annotation"]
    output:
        isoforms = config["output"]["transcript_dir"] + "/{sample}_isoforms.csv",
        counts = config["output"]["transcript_dir"] + "/{sample}_counts.csv"
    log:
        config["output"]["transcript_dir"] + "/logs/{sample}_transcript_id.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/transcript_identification.py"

rule compare_transcripts:
    input:
        tumor_isoforms = lambda wildcards: expand(
            config["output"]["transcript_dir"] + "/{sample}_isoforms.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "tumor"]
        ),
        tumor_counts = lambda wildcards: expand(
            config["output"]["transcript_dir"] + "/{sample}_counts.csv", 
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "tumor"]
        ),
        normal_isoforms = lambda wildcards: expand(
            config["output"]["transcript_dir"] + "/{sample}_isoforms.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "normal"]
        ),
        normal_counts = lambda wildcards: expand(
            config["output"]["transcript_dir"] + "/{sample}_counts.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "normal"]
        )
    output:
        differential = config["output"]["transcript_dir"] + "/differential_isoforms.csv",
        novel = config["output"]["transcript_dir"] + "/novel_isoforms.csv",
        summary = config["output"]["transcript_dir"] + "/isoform_summary.txt"
    log:
        config["output"]["transcript_dir"] + "/logs/transcript_comparison.log"
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/transcript_comparison.R"
