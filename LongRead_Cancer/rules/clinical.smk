rule integrate_clinical_data:
    input:
        isoforms = expand(
            config["output"]["transcript_dir"] + "/{sample}_isoforms.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "tumor"]
        ),
        junctions = expand(
            config["output"]["splicing_dir"] + "/{sample}_junctions.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "tumor"]
        ),
        fusions = expand(
            config["output"]["fusion_dir"] + "/{sample}_fusions.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "tumor"]
        ),
        ase = expand(
            config["output"]["ase_dir"] + "/{sample}_ase_results.csv",
            sample=[s for s in SAMPLES if config["samples"][s]["type"] == "tumor"]
        ),
        clinical = config["clinical_data"]
    output:
        integrated = config["output"]["clinical_dir"] + "/integrated_data.csv",
        survival = config["output"]["clinical_dir"] + "/survival_analysis.csv",
        plots_dir = directory(config["output"]["clinical_dir"] + "/plots")
    params:
        outcome_col = config["params"]["clinical"]["outcome_column"],
        event_col = config["params"]["clinical"]["event_column"],
        sample_id_col = config["params"]["clinical"]["sample_id_column"]
    log:
        config["output"]["clinical_dir"] + "/logs/clinical_integration.log"
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/clinical_integration.R"
