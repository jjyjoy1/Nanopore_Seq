# Nanopore Direct RNA Sequencing Pipeline for rRNA Modification Analysis

This Snakemake pipeline analyzes direct RNA sequencing data from Oxford Nanopore technology to identify ribosomal RNA (rRNA) modifications and classify tissue origin and tumor status based on modification patterns.

## Overview
Direct RNA sequencing allows the detection of RNA modifications by analyzing current disruptions as RNA molecules pass through nanopores. This pipeline processes raw FAST5 files from the sequencer, identifies rRNA modifications, and uses machine learning to classify samples.

## Prerequisites

Before running this pipeline, ensure you have the following installed:

- Python 3.7+
- Snakemake
- Oxford Nanopore software:
  - Guppy (for basecalling)
- Third-party tools:
  - minimap2 (for alignment)
  - samtools
  - Tombo (for modification detection)
- Python packages:
  - pandas
  - numpy
  - scikit-learn
  - matplotlib
  - seaborn

## Directory Structure

Set up your directory with the following structure:

```
project_directory/
├── Snakefile
├── config.yaml
├── scripts/
│   ├── prepare_ml_data.py
│   └── train_classifier.py
├── fast5_files/
│   ├── sample1/
│   ├── sample2/
│   └── ...
├── references/
│   └── transcriptome.fa
└── sample_metadata.csv
```

## Configuration

Edit the `config.yaml` file to specify:

1. Sample names
2. Reference files
3. rRNA IDs of interest
4. Directory paths
5. Computational resources

## Sample Metadata

The `sample_metadata.csv` file should contain at least the following columns:
- sample_id: Matches the sample directory names
- tissue_type: Type of tissue (e.g., lung, liver)
- cancer_status: "normal" or "cancer"

## Running the Pipeline

1. Install all prerequisites
2. Set up the directory structure
3. Configure the pipeline in `config.yaml`
4. Run the pipeline:

```bash
snakemake --cores 8
```

For a dry run (to check the execution plan):

```bash
snakemake --dry-run
```

To generate a workflow diagram:

```bash
snakemake --dag | dot -Tpng > dag.png
```

## Pipeline Steps

1. **Basecalling**: Convert raw signal data to nucleotide sequences
2. **Alignment**: Map reads to reference transcriptome
3. **Resquiggle**: Realign raw signal to reference sequence
4. **Modification Detection**: Identify modified nucleotides
5. **Data Preparation**: Process modification data for machine learning
6. **Classification**: Train model to classify tissue origin and tumor status

## Output

The pipeline generates:
- Processed modification data in CSV format
- A trained classifier model
- Classification report with performance metrics
- Visualization of important modification sites
- ROC curve and confusion matrix

## Customization

- To add custom analysis steps, modify the Snakefile or add new rules
- For different classification models, edit the `train_classifier.py` script
- To analyze different RNA types, change the reference transcriptome and target regions

## Troubleshooting

If you encounter issues:
1. Check log files in the output directory
2. Verify input file formats
3. Ensure all tools are properly installed
4. For memory issues, adjust threads in config.yaml

## References

For more information on nanopore direct RNA sequencing and modification analysis:
- Oxford Nanopore Technologies: [https://nanoporetech.com](https://nanoporetech.com)
- Tombo documentation: [https://nanoporetech.github.io/tombo/](https://nanoporetech.github.io/tombo/)

