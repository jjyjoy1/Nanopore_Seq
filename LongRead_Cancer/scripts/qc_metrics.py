#!/usr/bin/env python
"""
Calculate quality control metrics for long-read RNA-seq data
"""

import sys
import pandas as pd
import numpy as np
import pysam

def calculate_qc_metrics(bam_file, output_file, min_read_length=500):
    """
    Calculate basic QC metrics for long-read data
    """
    
    metrics = {}
    total_reads = 0
    mapped_reads = 0
    read_lengths = []
    
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Iterate through reads
    for read in bam:
        total_reads += 1
        
        # Store read length
        read_length = len(read.query_sequence) if read.query_sequence else 0
        if read_length >= min_read_length:
            read_lengths.append(read_length)
        
        # Check if read is mapped
        if not read.is_unmapped:
            mapped_reads += 1
    
    # Calculate mapping rate
    mapping_rate = (mapped_reads / total_reads) * 100 if total_reads > 0 else 0
    
    # Calculate read length statistics
    metrics['total_reads'] = total_reads
    metrics['mapped_reads'] = mapped_reads
    metrics['mapping_rate'] = mapping_rate
    metrics['median_read_length'] = np.median(read_lengths) if read_lengths else 0
    metrics['mean_read_length'] = np.mean(read_lengths) if read_lengths else 0
    metrics['min_read_length'] = np.min(read_lengths) if read_lengths else 0
    metrics['max_read_length'] = np.max(read_lengths) if read_lengths else 0
    metrics['reads_>10kb'] = sum(1 for l in read_lengths if l > 10000)
    metrics['reads_>5kb'] = sum(1 for l in read_lengths if l > 5000)
    
    # Write metrics to output file
    with open(output_file, 'w') as f:
        for key, value in metrics.items():
            f.write(f"{key}\t{value}\n")
    
    return metrics

# Snakemake integration
if __name__ == "__main__":
    bam_file = snakemake.input.bam
    output_file = snakemake.output.metrics
    min_read_length = snakemake.params.get("min_read_length", 500)
    
    calculate_qc_metrics(bam_file, output_file, min_read_length)


