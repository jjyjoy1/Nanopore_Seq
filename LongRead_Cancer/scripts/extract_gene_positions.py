#!/usr/bin/env python
"""
Extract gene positions from GTF file for fusion detection
"""

import sys
import pickle
import re
from collections import defaultdict
from intervaltree import IntervalTree

def extract_gene_positions(gtf_file, output_file):
    """
    Extract gene positions from GTF file and create interval trees
    """
    
    gene_positions = defaultdict(list)
    gene_names = {}
    
    # Parse GTF file
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if fields[2] != 'gene':
                continue
                
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            
            # Extract gene ID and name from attributes
            attributes = fields[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
            
            if gene_id_match:
                gene_id = gene_id_match.group(1)
                
                # Get gene name if available, otherwise use gene_id
                gene_name = gene_id
                if gene_name_match:
                    gene_name = gene_name_match.group(1)
                    
                gene_positions[chrom].append((start, end, gene_id))
                gene_names[gene_id] = gene_name
    
    # Create interval trees for fast overlap checking
    gene_trees = {}
    for chrom, positions in gene_positions.items():
        tree = IntervalTree()
        for start, end, gene_id in positions:
            tree[start:end] = gene_id
        gene_trees[chrom] = tree
    
    # Save to output file
    with open(output_file, 'wb') as f:
        pickle.dump((gene_trees, gene_names), f)
    
    return gene_trees, gene_names

# Snakemake integration
if __name__ == "__main__":
    gtf_file = snakemake.input.annotation
    output_file = snakemake.output.gene_positions
    
    extract_gene_positions(gtf_file, output_file)

