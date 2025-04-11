#!/usr/bin/env python
"""
Filter and analyze tumor-specific fusion candidates
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def filter_tumor_specific_fusions(tumor_fusions_files, normal_fusions_files, 
                                 output_file, plots_dir):
    """
    Filter for tumor-specific fusion candidates
    """
    
    # Create plots directory
    os.makedirs(plots_dir, exist_ok=True)
    
    # Combine tumor and normal fusion data
    tumor_fusions_all = []
    for file in tumor_fusions_files:
        fusions = pd.read_csv(file)
        sample_id = os.path.basename(file).split('_')[0]
        fusions['sample'] = sample_id
        tumor_fusions_all.append(fusions)
    
    normal_fusions_all = []
    for file in normal_fusions_files:
        fusions = pd.read_csv(file)
        sample_id = os.path.basename(file).split('_')[0]
        fusions['sample'] = sample_id
        normal_fusions_all.append(fusions)
    
    # Combine all samples
    if tumor_fusions_all:
        tumor_fusions = pd.concat(tumor_fusions_all)
    else:
        tumor_fusions = pd.DataFrame(columns=['gene1_id', 'gene2_id', 'supporting_reads',
                                             'gene1_name', 'gene2_name', 'fusion_name', 'sample'])
    
    if normal_fusions_all:
        normal_fusions = pd.concat(normal_fusions_all)
    else:
        normal_fusions = pd.DataFrame(columns=['gene1_id', 'gene2_id', 'supporting_reads',
                                             'gene1_name', 'gene2_name', 'fusion_name', 'sample'])
    
    # Create sets of fusion pairs in normal samples
    normal_fusion_pairs = set(zip(normal_fusions['gene1_id'], normal_fusions['gene2_id']))
    
    # Filter for tumor-specific fusions
    tumor_specific_fusions = tumor_fusions[~tumor_fusions.apply(
        lambda x: (x['gene1_id'], x['gene2_id']) in normal_fusion_pairs, axis=1
    )]
    
    # Aggregate reads across samples for each fusion
    if len(tumor_specific_fusions) > 0:
        tumor_specific_aggregated = tumor_specific_fusions.groupby(
            ['gene1_id', 'gene2_id', 'gene1_name', 'gene2_name', 'fusion_name']
        )['supporting_reads'].sum().reset_index()
    else:
        tumor_specific_aggregated = pd.DataFrame(columns=['gene1_id', 'gene2_id', 'supporting_reads',
                                                        'gene1_name', 'gene2_name', 'fusion_name'])
    
    # Save results
    tumor_specific_aggregated.to_csv(output_file, index=False)
    
    # Create visualizations
    
    # 1. Top fusions by supporting read count
    if len(tumor_specific_aggregated) > 0:
        top_fusions = tumor_specific_aggregated.sort_values('supporting_reads', ascending=False).head(20)
        
        plt.figure(figsize=(12, 8))
        plt.barh(range(len(top_fusions)), top_fusions['supporting_reads'], color='teal')
        plt.yticks(range(len(top_fusions)), top_fusions['fusion_name'])
        plt.xlabel('Supporting Reads')
        plt.title('Top 20 Tumor-Specific Fusion Candidates')
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/top_tumor_specific_fusions.png")
        plt.close()
        
        # 2. Distribution of supporting read counts
        plt.figure(figsize=(10, 6))
        plt.hist(tumor_specific_aggregated['supporting_reads'], bins=20, color='skyblue', edgecolor='black')
        plt.xlabel('Supporting Reads')
        plt.ylabel('Number of Fusion Candidates')
        plt.title('Distribution of Supporting Reads for Tumor-Specific Fusions')
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/fusion_support_distribution.png")
        plt.close()
        
        # 3. Network visualization of gene fusion relationships
        if len(tumor_specific_aggregated) > 0 and len(tumor_specific_aggregated) <= 50:
            plt.figure(figsize=(14, 14))
            
            # Create a network graph
            from matplotlib.lines import Line2D
            
            # Extract genes
            all_genes = set()
            gene_pairs = []
            
            for _, row in tumor_specific_aggregated.iterrows():
                all_genes.add(row['gene1_name'])
                all_genes.add(row['gene2_name'])
                gene_pairs.append((row['gene1_name'], row['gene2_name'], row['supporting_reads']))
            
            # Position genes in a circle
            import math
            
            n_genes = len(all_genes)
            gene_positions = {}
            
            for i, gene in enumerate(sorted(all_genes)):
                angle = 2 * math.pi * i / n_genes
                x = math.cos(angle)
                y = math.sin(angle)
                gene_positions[gene] = (x, y)
            
            # Draw genes as nodes
            for gene, (x, y) in gene_positions.items():
                plt.scatter(x, y, s=300, c='skyblue', edgecolor='black', zorder=10)
                plt.text(x * 1.1, y * 1.1, gene, fontsize=8, ha='center', va='center')
            
            # Draw connections between fused genes
            max_support = max([s for _, _, s in gene_pairs])
            
            for gene1, gene2, support in gene_pairs:
                x1, y1 = gene_positions[gene1]
                x2, y2 = gene_positions[gene2]
                
                # Scale line width by supporting reads
                line_width = 0.5 + 3 * (support / max_support)
                
                plt.plot([x1, x2], [y1, y2], 'r-', alpha=0.6, linewidth=line_width)
            
            # Add legend
            legend_elements = [
                Line2D([0], [0], color='r', linewidth=1, alpha=0.6, label='Low Support'),
                Line2D([0], [0], color='r', linewidth=2, alpha=0.6, label='Medium Support'),
                Line2D([0], [0], color='r', linewidth=3, alpha=0.6, label='High Support')
            ]
            plt.legend(handles=legend_elements, loc='upper right')
            
            # Remove axes
            plt.axis('off')
            plt.title('Gene Fusion Network')
            plt.tight_layout()
            plt.savefig(f"{plots_dir}/fusion_network.png")
            plt.close()
    
    return tumor_specific_aggregated

# Snakemake integration
if __name__ == "__main__":
    tumor_fusions_files = snakemake.input.tumor_fusions
    normal_fusions_files = snakemake.input.normal_fusions
    output_file = snakemake.output.tumor_specific
    plots_dir = snakemake.output.plots_dir
    
    filter_tumor_specific_fusions(tumor_fusions_files, normal_fusions_files, output_file, plots_dir)

