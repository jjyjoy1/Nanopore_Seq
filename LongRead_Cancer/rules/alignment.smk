rule align_reads:
    input:
        fastq = lambda wildcards: config["samples"][wildcards.sample]["fastq"],
        reference = config["reference"]["genome"]
    output:
        bam = config["output"]["align_dir"] + "/{sample}.bam",
        bai = config["output"]["align_dir"] + "/{sample}.bam.bai"
    log:
        config["output"]["align_dir"] + "/logs/{sample}_alignment.log"
    conda:
        "../envs/base.yaml"
    threads: 
        config["params"]["alignment"]["threads"]
    shell:
        """
        # Align with minimap2 and pipe to samtools for sorting
        minimap2 -ax splice:hq {config[params][alignment][extra]} \
            -t {threads} {input.reference} {input.fastq} | \
            samtools sort -@ {threads} -o {output.bam} \
            2> {log}
            
        # Index the BAM file
        samtools index {output.bam}
        """
