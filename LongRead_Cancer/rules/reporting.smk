rule create_qc_report:
    input:
        qc_metrics = expand(config["output"]["qc_dir"] + "/{sample}_alignment_metrics.txt", sample=SAMPLES),
        nanoplot = expand(config["output"]["qc_dir"] + "/{sample}_nanoplot.html", sample=SAMPLES),
        multiqc = config["output"]["qc_dir"] + "/multiqc_report.html"
    output:
        report = config["output"]["qc_dir"] + "/qc_report.html"
    log:
        config["output"]["qc_dir"] + "/logs/create_qc_report.log"
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/create_qc_report.Rmd"

rule create_transcript_report:
    input:
        differential = config["output"]["transcript_dir"] + "/differential_isoforms.csv",
        novel = config["output"]["transcript_dir"] + "/novel_isoforms.csv",
        summary = config["output"]["transcript_dir"] + "/isoform_summary.txt"
    output:
        report = config["output"]["transcript_dir"] + "/transcript_summary.html"
    log:
        config["output"]["transcript_dir"] + "/logs/create_transcript_report.log"
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/create_transcript_report.Rmd"

rule create_splicing_report:
    input:
        tumor_specific = config["output"]["splicing_dir"] + "/tumor_specific_junctions.csv",
        differential = config["output"]["splicing_dir"] + "/differential_junctions.csv",
        plots_dir = config["output"]["splicing_dir"] + "/plots"
    output:
        report = config["output"]["splicing_dir"] + "/splicing_analysis.html"
    log:
        config["output"]["splicing_dir"] + "/logs/create_splicing_report.log"
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/create_splicing_report.Rmd"

rule create_fusion_report:
    input:
        tumor_specific = config["output"]["fusion_dir"] + "/tumor_specific_fusions.csv",
        plots_dir = config["output"]["fusion_dir"] + "/plots"
    output:
        report = config["output"]["fusion_dir"] + "/fusion_detection.html"
    log:
        config["output"]["fusion_dir"] + "/logs/create_fusion_report.log"
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/create_fusion_report.Rmd"

rule create_ase_report:
    input:
        ase_results = expand(
            config["output"]["ase_dir"] + "/{sample}_ase_results.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "tumor"]
        ),
        plots_dir = expand(
            config["output"]["ase_dir"] + "/{sample}_ase_plots",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "tumor"]
        )
    output:
        report = config["output"]["ase_dir"] + "/ase_analysis.html"
    log:
        config["output"]["ase_dir"] + "/logs/create_ase_report.log"
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/create_ase_report.Rmd"

rule create_clinical_report:
    input:
        integrated = config["output"]["clinical_dir"] + "/integrated_data.csv",
        survival = config["output"]["clinical_dir"] + "/survival_analysis.csv",
        plots_dir = config["output"]["clinical_dir"] + "/plots"
    output:
        report = config["output"]["clinical_dir"] + "/clinical_integration.html"
    log:
        config["output"]["clinical_dir"] + "/logs/create_clinical_report.log"
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/create_clinical_report.Rmd"


