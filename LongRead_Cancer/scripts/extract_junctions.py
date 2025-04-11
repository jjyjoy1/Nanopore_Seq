#!/usr/bin/env python
"""
Extract splice junctions from long-read alignments
"""

import sys
import pandas as pd
from collections import defaultdict
import pysam

def extract_splice_junctions(bam_file, output_file, min_support=3):
    """
    Extract splice junctions from a BAM file
    """
    
    junctions = defaultdict(int)
    
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Iterate through aligned reads
    for read in bam:
        if read.is_unmapped:
            continue
            
        # Skip reads without CIGAR string
        if not read.cigartuples:
            continue
            
        # Extract positions and reference positions
        pos = read.pos
        chrom = bam.get_reference_name(read.reference_id)
        
        # Analyze CIGAR string to find splice junctions (N operations)
        for operation, length in read.cigartuples:
            if operation == 0:  # Match
                pos += length
            elif operation == 1:  # Insertion
                continue
            elif operation == 2:  # Deletion
                pos += length
            elif operation == 3:  # Splice junction (N)
                junction_start = pos
                junction_end = pos + length
                junction_id = f"{chrom}:{junction_start}-{junction_end}"
                junctions[junction_id] += 1
                pos += length
                
    # Filter junctions by minimum support
    filtered_junctions = {j: count for j, count in junctions.items() if count >= min_support}
    
    # Convert to DataFrame
    junctions_df = pd.DataFrame({
        'junction_id': list(filtered_junctions.keys()),
        'read_count': list(filtered_junctions.values())
    })
    
    # Parse junction coordinates
    junctions_df[['chrom', 'coordinates']] = junctions_df['junction_id'].str.split(':', expand=True)
    junctions_df[['start', 'end']] = junctions_df['coordinates'].str.split('-', expand=True)
    junctions_df['start'] = junctions_df['start'].astype(int)
    junctions_df['end'] = junctions_df['end'].astype(int)
    
    # Write to output file
    junctions_df.to_csv(output_file, index=False)
    
    return junctions_df

# Snakemake integration
if __name__ == "__main__":
    bam_file = snakemake.input.bam
    output_file = snakemake.output.junctions
    min_support = snakemake.params.min_support
    
    extract_splice_junctions(bam_file, output_file, min_support)

