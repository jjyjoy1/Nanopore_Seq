#!/usr/bin/env python
"""
Identify and quantify full-length transcripts from long reads
"""

import os
import sys
import json
import subprocess
import pandas as pd

def identify_transcripts(bam_file, reference_genome, reference_gtf, 
                        output_isoforms, output_counts, min_coverage=3):
    """
    Identify and quantify full-length transcripts using FLAMES
    """
    
    # Create temporary output directory
    output_dir = os.path.dirname(output_isoforms)
    os.makedirs(output_dir, exist_ok=True)
    
    # Get sample name from bam file
    sample_name = os.path.basename(bam_file).split('.')[0]
    
    # Set up FLAMES config
    config_file = f"{output_dir}/{sample_name}_config.json"
    
    config = {
        "genome_annotation": reference_gtf,
        "genome_dir": os.path.dirname(reference_genome),
        "genome_file": os.path.basename(reference_genome),
        "outdir": output_dir,
        "thread": 16,
        "sample_name": sample_name,
        "fastq_file": "",  # We're using the BAM directly
        "bam_file": bam_file,
        "do_genome_alignment": False,  # Already aligned
        "minimap2_dir": "",  # Will use system installation
        "annotation_flatten_flag": True,
        "transcript_counting": True,
        "do_isoform_identification": True,
        "do_read_realignment": True,
        "do_error_correction": True,
        "min_transcript_coverage": min_coverage
    }
    
    # Write config to file
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=4)
    
    # Run FLAMES
    subprocess.run(f"flames {config_file}", shell=True, check=True)
    
    # Copy results to output files (FLAMES creates its own naming scheme)
    flames_isoforms = f"{output_dir}/{sample_name}_isoforms.csv"
    flames_counts = f"{output_dir}/{sample_name}_counts.csv"
    
    # If files exist, copy them to the specified output
    if os.path.exists(flames_isoforms):
        subprocess.run(f"cp {flames_isoforms} {output_isoforms}", shell=True, check=True)
    else:
        # Create empty files if FLAMES didn't produce output
        pd.DataFrame().to_csv(output_isoforms, index=False)
        
    if os.path.exists(flames_counts):
        subprocess.run(f"cp {flames_counts} {output_counts}", shell=True, check=True)
    else:
        pd.DataFrame().to_csv(output_counts, index=False)

# Snakemake integration
if __name__ == "__main__":
    bam_file = snakemake.input.bam
    reference_genome = snakemake.input.reference
    reference_gtf = snakemake.input.annotation
    output_isoforms = snakemake.output.isoforms
    output_counts = snakemake.output.counts
    min_coverage = snakemake.params.get("min_coverage", 3)
    
    identify_transcripts(bam_file, reference_genome, reference_gtf, 
                         output_isoforms, output_counts, min_coverage)

