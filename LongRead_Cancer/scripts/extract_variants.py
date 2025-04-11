#!/usr/bin/env python
"""
Extract variants from long reads for allele-specific expression analysis
"""

import sys
import pickle
from collections import defaultdict
import pysam

def extract_variants_from_reads(bam_file, vcf_file, output_file, region=None):
    """
    Extract variants from long reads and store read-variant associations
    """
    
    # Open BAM and VCF files
    bam = pysam.AlignmentFile(bam_file, "rb")
    vcf = pysam.VariantFile(vcf_file, "r")
    
    # Dictionary to store read-variant associations
    read_variants = defaultdict(list)
    variant_positions = {}
    
    # Process VCF to get variant positions
    for variant in vcf.fetch(region=region):
        var_id = f"{variant.chrom}:{variant.pos}_{variant.ref}_{','.join(variant.alts)}"
        variant_positions[(variant.chrom, variant.pos)] = var_id
    
    # Process reads to find variants
    for read in bam.fetch(region=region):
        if read.is_unmapped:
            continue
            
        read_id = read.query_name
        read_seq = read.query_sequence
        
        # Get reference positions
        aligned_pairs = read.get_aligned_pairs(with_seq=True)
        
        for read_pos, ref_pos, ref_base in aligned_pairs:
            # Skip positions without alignment
            if read_pos is None or ref_pos is None:
                continue
                
            chrom = bam.get_reference_name(read.reference_id)
            var_key = (chrom, ref_pos + 1)  # Convert to 1-based for VCF compatibility
            
            # Check if position is a known variant
            if var_key in variant_positions:
                read_base = read_seq[read_pos]
                var_id = variant_positions[var_key]
                
                # Store variant info for this read
                read_variants[read_id].append((var_id, read_base))
    
    # Save to output file
    with open(output_file, 'wb') as f:
        pickle.dump(read_variants, f)
    
    return read_variants

# Snakemake integration
if __name__ == "__main__":
    bam_file = snakemake.input.bam
    vcf_file = snakemake.input.vcf
    output_file = snakemake.output.variants
    region = snakemake.params.get("region", None)
    
    extract_variants_from_reads(bam_file, vcf_file, output_file, region)

