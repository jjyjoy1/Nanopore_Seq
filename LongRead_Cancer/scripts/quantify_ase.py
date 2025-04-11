#!/usr/bin/env python
"""
Quantify allele-specific expression using phased variants
"""

import sys
import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pysam
import re

def extract_gene_regions(gtf_file, genes_of_interest=None):
    """
    Extract genomic regions for genes of interest
    """
    
    gene_regions = {}
    
    # Parse GTF file
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if fields[2] != 'gene':
                continue
                
            # Extract attributes
            attributes = fields[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            
            if not gene_id_match:
                continue
                
            gene_id = gene_id_match.group(1)
            
            # Skip if not in genes of interest (if specified)
            if genes_of_interest and gene_id not in genes_of_interest:
                continue
                
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            
            gene_regions[gene_id] = (chrom, start, end)
    
    return gene_regions

def quantify_ase(bam_file, phased_file, gtf_file, output_file, plots_dir, min_coverage=10,
                genes_of_interest=None):
    """
    Quantify allele-specific expression using phased variants
    """
    
    # Create plots directory
    os.makedirs(plots_dir, exist_ok=True)
    
    # Load phased blocks
    with open(phased_file, 'rb') as f:
        phased_blocks = pickle.load(f)
    
    # Extract gene regions
    gene_regions = extract_gene_regions(gtf_file, genes_of_interest)
    
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    ase_results = []
    
    # Process each gene region
    for gene_id, region in gene_regions.items():
        chrom, start, end = region
        
        # Count reads supporting each haplotype
        hap1_count = 0
        hap2_count = 0
        unphased_count = 0
        
        # Process reads in this region
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped:
                continue
                
            read_id = read.query_name
            
            # Try to assign this read to a haplotype
            assigned = False
            
            # Check variants in this read against phased blocks
            for block in phased_blocks:
                hap1_matches = 0
                hap2_matches = 0
                
                # Extract alleles from read sequence
                read_seq = read.query_sequence
                aligned_pairs = read.get_aligned_pairs(with_seq=True)
                read_variants = {}
                
                for read_pos, ref_pos, ref_base in aligned_pairs:
                    if read_pos is None or ref_pos is None:
                        continue
                        
                    var_key = f"{chrom}:{ref_pos+1}"
                    read_allele = read_seq[read_pos]
                    read_variants[var_key] = read_allele
                
                # Check if read variants match phased variants
                for var_id, alleles in block['variants'].items():
                    var_prefix = var_id.split('_')[0]  # Extract position part
                    
                    if var_prefix in read_variants:
                        read_allele = read_variants[var_prefix]
                        
                        if alleles[0] == read_allele:
                            hap1_matches += 1
                        elif len(alleles) > 1 and alleles[1] == read_allele:
                            hap2_matches += 1
                
                # Assign to haplotype with more matches
                if hap1_matches > 0 or hap2_matches > 0:
                    if hap1_matches > hap2_matches:
                        hap1_count += 1
                        assigned = True
                        break
                    elif hap2_matches > hap1_matches:
                        hap2_count += 1
                        assigned = True
                        break
            
            # Count unphased reads
            if not assigned:
                unphased_count += 1
        
        # Calculate allelic ratio if coverage is sufficient
        total_phased = hap1_count + hap2_count
        if total_phased >= min_coverage:
            hap1_ratio = hap1_count / total_phased
            hap2_ratio = hap2_count / total_phased
            ase_score = abs(hap1_ratio - 0.5) * 2  # 0=balanced, 1=monoallelic
        else:
            hap1_ratio = np.nan
            hap2_ratio = np.nan
            ase_score = np.nan
        
        # Store results
        ase_results.append({
            'gene_id': gene_id,
            'haplotype1_count': hap1_count,
            'haplotype2_count': hap2_count,
            'unphased_count': unphased_count,
            'total_phased': total_phased,
            'haplotype1_ratio': hap1_ratio,
            'haplotype2_ratio': hap2_ratio,
            'ase_score': ase_score
        })
        
        # Create visualization for this gene
        if total_phased >= min_coverage:
            plt.figure(figsize=(8, 6))
            plt.bar(['Haplotype 1', 'Haplotype 2'], [hap1_count, hap2_count], color=['blue', 'orange'])
            plt.title(f'Allele-Specific Expression for {gene_id}')
            plt.ylabel('Read Count')
            plt.savefig(f"{plots_dir}/{gene_id}_ase.png")
            plt.close()
    
    # Convert to DataFrame and save
    ase_df = pd.DataFrame(ase_results)
    ase_df.to_csv(output_file, index=False)
    
    # Create summary plot of ASE scores
    if len(ase_df) > 0 and not ase_df['ase_score'].isna().all():
        plt.figure(figsize=(10, 6))
        sorted_df = ase_df.sort_values('ase_score', ascending=False).dropna(subset=['ase_score'])
        plt.bar(sorted_df['gene_id'], sorted_df['ase_score'], color='purple')
        plt.axhline(y=0.3, color='red', linestyle='--', alpha=0.7)  # Typical threshold for ASE
        plt.xticks(rotation=90)
        plt.title('Allele-Specific Expression Scores')
        plt.ylabel('ASE Score (0 = balanced, 1 = monoallelic)')
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/ase_scores_summary.png")
        plt.close()
    
    return ase_df

# Snakemake integration
if __name__ == "__main__":
    bam_file = snakemake.input.bam
    phased_file = snakemake.input.phased
    gtf_file = snakemake.input.annotation
    output_file = snakemake.output.ase_results
    plots_dir = snakemake.output.plots_dir
    min_coverage = snakemake.params.min_coverage
    
    # Define a list of cancer-related genes (optional)
    cancer_genes = [
        'TP53', 'BRCA1', 'BRCA2', 'EGFR', 'KRAS', 'PIK3CA', 'PTEN', 'RB1', 'CDKN2A',
        'APC', 'BRAF', 'NRAS', 'IDH1', 'IDH2', 'KIT', 'PDGFRA', 'FGFR1', 'FGFR2',
        'FGFR3', 'ALK', 'ROS1', 'RET', 'ERBB2', 'MYC', 'MYCN', 'CDH1', 'CTNNB1'
    ]
    
    quantify_ase(bam_file, phased_file, gtf_file, output_file, plots_dir, min_coverage, cancer_genes)

