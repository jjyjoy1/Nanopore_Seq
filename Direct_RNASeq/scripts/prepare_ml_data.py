#!/usr/bin/env python
# prepare_ml_data.py - Process modification data for ML
# Save in scripts/prepare_ml_data.py

import pandas as pd
import os
import numpy as np
from snakemake.shell import shell

# Get input and output file paths from Snakemake
mod_file = snakemake.input.mods
metadata_file = snakemake.input.metadata
output_file = snakemake.output.ml_data

# Read modification data
print(f"Reading modification data from {mod_file}")
mod_data = pd.read_csv(mod_file, sep='\t')

# Read sample metadata
print(f"Reading sample metadata from {metadata_file}")
metadata = pd.read_csv(metadata_file)

# Process modification data
# The format from Tombo is typically: ref_id, pos, fraction_modified, ...
# We need to restructure this to have samples as rows and modification sites as columns

# Get unique sample IDs from the metadata
sample_ids = metadata['sample_id'].unique()

# For each sample, we need to extract its modification data
# This assumes Tombo outputs contain sample information
# If not, we'll need to process each sample separately and combine

# Check if mod_data has a sample column
if 'sample_id' not in mod_data.columns:
    print("Modification data doesn't contain sample_id column.")
    print("Processing each sample's modification files separately and combining...")
    
    # Get the directory containing the mod_file
    mod_dir = os.path.dirname(mod_file)
    
    # Create an empty list to store dataframes for each sample
    sample_dfs = []
    
    # Process each sample's modification file
    for sample_id in sample_ids:
        # Assuming each sample has its own file with naming pattern: sample_id.dampened_fraction.csv
        sample_mod_file = os.path.join(mod_dir, f"{sample_id}.dampened_fraction.csv")
        
        if os.path.exists(sample_mod_file):
            # Read the sample's modification data
            sample_mod_data = pd.read_csv(sample_mod_file, sep='\t')
            
            # Add sample_id column
            sample_mod_data['sample_id'] = sample_id
            
            # Append to the list
            sample_dfs.append(sample_mod_data)
        else:
            print(f"Warning: No modification file found for sample {sample_id}")
    
    # Combine all sample dataframes
    if sample_dfs:
        mod_data = pd.concat(sample_dfs, ignore_index=True)
    else:
        raise ValueError("No sample modification data could be processed")
else:
    print("Using sample_id column from modification data")

# Create a pivot table: rows=samples, columns=modification sites
print("Creating feature matrix of modification sites...")
mod_matrix = mod_data.pivot_table(
    index='sample_id',
    columns=['ref_id', 'pos'],
    values='fraction_modified',
    fill_value=0
)

# Flatten the multi-index columns
mod_matrix.columns = [f"{ref_id}_{pos}" for ref_id, pos in mod_matrix.columns]

# Merge with metadata
print("Merging with metadata...")
full_dataset = pd.merge(metadata, mod_matrix, on='sample_id', how='inner')

# Check if any samples were lost in the merge
merged_samples = full_dataset['sample_id'].unique()
missing_samples = set(sample_ids) - set(merged_samples)
if missing_samples:
    print(f"Warning: The following samples were not found in the merged data: {missing_samples}")

# Save the final dataset
print(f"Saving processed data to {output_file}")
full_dataset.to_csv(output_file, index=False)

print(f"Created ML dataset with {full_dataset.shape[0]} samples and {full_dataset.shape[1]} features")
print(f"Feature columns include: {', '.join(full_dataset.columns[:10])}...")

