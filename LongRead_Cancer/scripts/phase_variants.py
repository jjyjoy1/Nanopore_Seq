#!/usr/bin/env python
"""
Phase variants based on co-occurrence in long reads
"""

import sys
import pickle
from collections import defaultdict, deque

def phase_variants(variants_file, output_file, min_support=2):
    """
    Group variants into haplotypes based on co-occurrence in reads
    """
    
    # Load variants
    with open(variants_file, 'rb') as f:
        read_variants = pickle.load(f)
    
    # Count co-occurrences of variant pairs
    pair_counts = defaultdict(int)
    
    for read_id, variants in read_variants.items():
        # Skip reads with only one variant
        if len(variants) < 2:
            continue
            
        # Count all pairs of variants in this read
        for i in range(len(variants)):
            for j in range(i+1, len(variants)):
                var1_id, var1_allele = variants[i]
                var2_id, var2_allele = variants[j]
                
                pair_key = ((var1_id, var1_allele), (var2_id, var2_allele))
                pair_counts[pair_key] += 1
    
    # Filter by minimum support
    filtered_pairs = {pair: count for pair, count in pair_counts.items() 
                     if count >= min_support}
    
    # Create a graph representation for connected components
    graph = defaultdict(set)
    for ((var1_id, var1_allele), (var2_id, var2_allele)), _ in filtered_pairs.items():
        graph[(var1_id, var1_allele)].add((var2_id, var2_allele))
        graph[(var2_id, var2_allele)].add((var1_id, var1_allele))
    
    # Find connected components (phased blocks)
    phased_blocks = []
    visited = set()
    
    for node in graph:
        if node in visited:
            continue
            
        # BFS to find all connected nodes
        block = set()
        queue = deque([node])
        
        while queue:
            current = queue.popleft()
            if current in visited:
                continue
                
            visited.add(current)
            block.add(current)
            
            for neighbor in graph[current]:
                if neighbor not in visited:
                    queue.append(neighbor)
        
        phased_blocks.append(block)
    
    # Convert blocks to a more readable format
    phased_results = []
    
    for block_idx, block in enumerate(phased_blocks):
        block_variants = defaultdict(list)
        
        for var_id, allele in block:
            block_variants[var_id].append(allele)
            
        phased_results.append({
            'block_id': block_idx,
            'variants': dict(block_variants),
            'size': len(block_variants)
        })
    
    # Save to output file
    with open(output_file, 'wb') as f:
        pickle.dump(phased_results, f)
    
    return phased_results

# Snakemake integration
if __name__ == "__main__":
    variants_file = snakemake.input.variants
    output_file = snakemake.output.phased
    min_support = snakemake.params.get("min_support", 2)
    
    phase_variants(variants_file, output_file, min_support)

