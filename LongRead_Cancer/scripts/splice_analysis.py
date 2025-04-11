#!/usr/bin/env python
"""
Compare splice junctions between tumor and normal samples
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def compare_junctions(tumor_junctions_files, normal_junctions_files, 
                      tumor_specific_output, differential_output, plots_dir):
    """
    Compare splice junctions between tumor and normal samples
    """
    
    # Create plots directory
    os.makedirs(plots_dir, exist_ok=True)
    
    # Combine tumor and normal junction data
    tumor_junctions_all = []
    for file in tumor_junctions_files:
        sample_junctions = pd.read_csv(file)
        sample_id = os.path.basename(file).split('_')[0]
        sample_junctions['sample'] = sample_id
        tumor_junctions_all.append(sample_junctions)
    
    normal_junctions_all = []
    for file in normal_junctions_files:
        sample_junctions = pd.read_csv(file)
        sample_id = os.path.basename(file).split('_')[0]
        sample_junctions['sample'] = sample_id
        normal_junctions_all.append(sample_junctions)
    
    # Combine all samples
    if tumor_junctions_all:
        tumor_junctions = pd.concat(tumor_junctions_all)
        # Get sum of read counts per junction across samples
        tumor_junctions_summed = tumor_junctions.groupby('junction_id')['read_count'].sum().reset_index()
    else:
        tumor_junctions_summed = pd.DataFrame(columns=['junction_id', 'read_count'])
    
    if normal_junctions_all:
        normal_junctions = pd.concat(normal_junctions_all)
        # Get sum of read counts per junction across samples
        normal_junctions_summed = normal_junctions.groupby('junction_id')['read_count'].sum().reset_index()
    else:
        normal_junctions_summed = pd.DataFrame(columns=['junction_id', 'read_count'])
    
    # Create sets of junction IDs
    tumor_junction_set = set(tumor_junctions_summed['junction_id'])
    normal_junction_set = set(normal_junctions_summed['junction_id'])
    
    # Find tumor-specific junctions
    tumor_specific = tumor_junction_set - normal_junction_set
    tumor_specific_df = tumor_junctions_summed[tumor_junctions_summed['junction_id'].isin(tumor_specific)]
    
    # Parse junction coordinates for tumor-specific junctions
    if len(tumor_specific_df) > 0 and 'chrom' not in tumor_specific_df.columns:
        tumor_specific_df[['chrom', 'coordinates']] = tumor_specific_df['junction_id'].str.split(':', expand=True)
        tumor_specific_df[['start', 'end']] = tumor_specific_df['coordinates'].str.split('-', expand=True)
    
    # Find shared junctions with differential usage
    shared_junctions = tumor_junction_set.intersection(normal_junction_set)
    
    # Create a comparison dataframe for shared junctions
    if shared_junctions:
        shared_comparison = pd.DataFrame({'junction_id': list(shared_junctions)})
        shared_comparison = shared_comparison.merge(
            tumor_junctions_summed[['junction_id', 'read_count']].rename(columns={'read_count': 'tumor_count'}),
            on='junction_id'
        )
        shared_comparison = shared_comparison.merge(
            normal_junctions_summed[['junction_id', 'read_count']].rename(columns={'read_count': 'normal_count'}),
            on='junction_id'
        )
        
        # Calculate fold change
        shared_comparison['fold_change'] = shared_comparison['tumor_count'] / shared_comparison['normal_count']
        shared_comparison['log2_fold_change'] = np.log2(shared_comparison['fold_change'])
        
        # Find differentially used junctions (simple threshold-based approach)
        differential_junctions = shared_comparison[abs(shared_comparison['log2_fold_change']) > 1]
    else:
        shared_comparison = pd.DataFrame(columns=['junction_id', 'tumor_count', 'normal_count', 
                                                 'fold_change', 'log2_fold_change'])
        differential_junctions = shared_comparison.copy()
    
    # Save results
    tumor_specific_df.to_csv(tumor_specific_output, index=False)
    differential_junctions.to_csv(differential_output, index=False)
    
    # Create visualizations
    
    # 1. Scatter plot comparing junction usage
    if len(shared_comparison) > 0:
        plt.figure(figsize=(10, 8))
        plt.scatter(shared_comparison['normal_count'], shared_comparison['tumor_count'], 
                   alpha=0.6, s=30, edgecolor='k', linewidth=0.5)
        
        # Add diagonal line
        max_val = max(shared_comparison['normal_count'].max(), shared_comparison['tumor_count'].max())
        plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.3)
        
        # Highlight differential junctions
        if len(differential_junctions) > 0:
            plt.scatter(differential_junctions['normal_count'], differential_junctions['tumor_count'],
                       alpha=0.8, s=50, color='red', edgecolor='k', linewidth=0.5)
        
        # Labels and title
        plt.xlabel('Normal Junction Count')
        plt.ylabel('Tumor Junction Count')
        plt.title('Splice Junction Usage: Tumor vs Normal')
        plt.xscale('log')
        plt.yscale('log')
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/junction_usage_comparison.png")
        plt.close()
        
        # 2. Histogram of fold changes
        plt.figure(figsize=(10, 6))
        plt.hist(shared_comparison['log2_fold_change'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        plt.axvline(x=0, color='red', linestyle='--')
        plt.xlabel('Log2 Fold Change (Tumor/Normal)')
        plt.ylabel('Number of Junctions')
        plt.title('Distribution of Splice Junction Usage Changes')
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/junction_fold_change_distribution.png")
        plt.close()
    
    # 3. Top 20 tumor-specific junctions
    if len(tumor_specific_df) > 0:
        top_tumor_specific = tumor_specific_df.sort_values('read_count', ascending=False).head(20)
        plt.figure(figsize=(12, 8))
        plt.barh(np.arange(len(top_tumor_specific)), top_tumor_specific['read_count'], color='firebrick')
        plt.yticks(np.arange(len(top_tumor_specific)), top_tumor_specific['junction_id'])
        plt.xlabel('Read Count')
        plt.title('Top 20 Tumor-Specific Splice Junctions')
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/top_tumor_specific_junctions.png")
        plt.close()
    
    return tumor_specific_df, differential_junctions

# Snakemake integration
if __name__ == "__main__":
    tumor_junctions_files = snakemake.input.tumor_junctions
    normal_junctions_files = snakemake.input.normal_junctions
    tumor_specific_output = snakemake.output.tumor_specific
    differential_output = snakemake.output.differential
    plots_dir = snakemake.output.plots_dir
    
    compare_junctions(tumor_junctions_files, normal_junctions_files, 
                     tumor_specific_output, differential_output, plots_dir)

