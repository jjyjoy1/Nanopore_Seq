#!/usr/bin/env python
"""
Detect fusion transcripts from long-read alignments
"""

import sys
import pickle
import pandas as pd
from collections import defaultdict
import pysam

def detect_fusion_candidates(bam_file, gene_positions_file, output_file, min_support=2):
    """
    Detect potential fusion transcripts from long reads
    """
    
    # Load gene positions (from extract_gene_positions.py)
    with open(gene_positions_file, 'rb') as f:
        gene_trees, gene_names = pickle.load(f)
    
    fusion_candidates = defaultdict(int)
    
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Process each read
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
            
        # Skip reads with low mapping quality
        if read.mapping_quality < 20:
            continue
            
        # Check for supplementary alignments (possible fusion indicator)
        if not read.has_tag('SA'):
            continue
            
        # Get primary alignment information
        primary_chrom = bam.get_reference_name(read.reference_id)
        primary_pos = read.reference_start
        
        # Get supplementary alignment information
        sa_tag = read.get_tag('SA')
        sa_alignments = sa_tag.strip().split(';')
        
        for sa_aln in sa_alignments:
            if not sa_aln:
                continue
                
            # Parse supplementary alignment
            sa_parts = sa_aln.split(',')
            if len(sa_parts) < 6:
                continue
                
            sa_chrom, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm = sa_parts
            sa_pos = int(sa_pos)
            
            # Skip same chromosome nearby alignments (likely structural misalignments, not fusions)
            if sa_chrom == primary_chrom and abs(sa_pos - primary_pos) < 100000:
                continue
                
            # Find genes that overlap with alignment positions
            primary_genes = []
            if primary_chrom in gene_trees:
                primary_overlaps = gene_trees[primary_chrom][primary_pos:primary_pos+1]
                primary_genes = [overlap.data for overlap in primary_overlaps]
            
            sa_genes = []
            if sa_chrom in gene_trees:
                sa_overlaps = gene_trees[sa_chrom][sa_pos:sa_pos+1]
                sa_genes = [overlap.data for overlap in sa_overlaps]
            
            # If both alignments overlap genes, potential fusion
            for gene1 in primary_genes:
                for gene2 in sa_genes:
                    if gene1 != gene2:  # Avoid self-fusions
                        # Create fusion pair (order alphabetically for consistency)
                        fusion = tuple(sorted([gene1, gene2]))
                        fusion_candidates[fusion] += 1
    
    # Filter by minimum support
    filtered_fusions = {fusion: count for fusion, count in fusion_candidates.items() 
                       if count >= min_support}
    
    # Convert to DataFrame
    fusions_df = pd.DataFrame({
        'gene1_id': [f[0] for f in filtered_fusions.keys()],
        'gene2_id': [f[1] for f in filtered_fusions.keys()],
        'supporting_reads': list(filtered_fusions.values())
    })
    
    # Add gene names
    fusions_df['gene1_name'] = fusions_df['gene1_id'].map(gene_names)
    fusions_df['gene2_name'] = fusions_df['gene2_id'].map(gene_names)
    
    # Fusion name
    fusions_df['fusion_name'] = fusions_df['gene1_name'] + '--' + fusions_df['gene2_name']
    
    # Write to output file
    fusions_df.to_csv(output_file, index=False)
    
    return fusions_df

# Snakemake integration
if __name__ == "__main__":
    bam_file = snakemake.input.bam
    gene_positions_file = snakemake.input.gene_positions
    output_file = snakemake.output.fusions
    min_support = snakemake.params.min_support
    
    detect_fusion_candidates(bam_file, gene_positions_file, output_file, min_support)


