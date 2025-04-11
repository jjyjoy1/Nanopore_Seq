# LongReadCancer: A Snakemake Pipeline for Long-Read RNA-seq in Cancer Genomics

A comprehensive Snakemake pipeline for analyzing long-read RNA sequencing data in cancer genomics, focusing on transcript isoform detection, alternative splicing, fusion transcript identification, and allele-specific expression.

## 🌟 Features

- **Full-length transcript analysis**: Identify and quantify complete transcript isoforms without assembly
- **Alternative splicing detection**: Discover cancer-specific splicing events with precise junction boundaries
- **Structural variant analysis**: Detect gene fusions and other complex rearrangements
- **Allele-specific expression**: Quantify haplotype-specific expression using long-read phasing
- **Clinical integration**: Correlate molecular features with patient outcomes
- **Comprehensive reporting**: Generate interactive HTML reports and visualizations
- **Reproducible**: Fully containerized workflow with versioned dependencies
- **Scalable**: Runs on workstations or HPC clusters with automatic parallelization

## 🔧 Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://github.com/mamba-org/mamba)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) ≥ 7.0.0

### Quick Start

1. Clone this repository:
```bash
git clone https://github.com/yourusername/LongReadCancer.git
cd LongReadCancer
```

2. Create the conda environments:
```bash
# Create all required environments at once
snakemake --use-conda --conda-create-envs-only
```

3. Edit the configuration file:
```bash
# Copy and modify the example config
cp config/config.example.yaml config/config.yaml
nano config/config.yaml
```

4. Run a test on the example data:
```bash
snakemake --use-conda -j 4 --configfile config/config.yaml
```

### Basic Usage

```bash
# Run the complete pipeline
snakemake --use-conda -j <cores> --configfile config/config.yaml

# Generate a workflow report
snakemake --report report.html

# Create a rulegraph visualization
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

### Running on a cluster

```bash
# Example for SLURM
snakemake --use-conda -j 100 --cluster "sbatch --mem={resources.mem_mb} -c {threads} -t {resources.runtime}"
```

### Module-specific analysis
```bash
# Run only the transcript identification
snakemake --use-conda -j <cores> transcript_report

# Run only the splicing analysis
snakemake --use-conda -j <cores> splicing_report

# Run only the fusion detection
snakemake --use-conda -j <cores> fusion_report

# Run only the ASE analysis
snakemake --use-conda -j <cores> ase_report
```

## 📂 Pipeline Components
### Input Data
- Long-read RNA-seq FASTQ files (ONT or PacBio)
- Reference genome FASTA
- Gene annotation GTF
- (Optional) Variant calls for ASE analysis
- (Optional) Clinical data for outcome analysis

### Analysis Modules

1. **Quality Control**
   - Raw read QC (NanoPlot)
   - Alignment quality metrics
   - MultiQC reports

2. **Transcript Identification**
   - Full-length transcript detection (FLAMES)
   - Novel isoform discovery
   - Differential isoform usage

3. **Splicing Analysis**
   - Splice junction extraction
   - Tumor-specific junction identification
   - Differential junction usage

4. **Fusion Transcript Detection**
   - Fusion candidate identification
   - Filtering against normal samples
   - Gene partner network analysis

5. **Allele-Specific Expression**
   - Variant extraction from long reads
   - Haplotype phasing
   - Allele-specific expression quantification

6. **Clinical Integration**
   - Feature-outcome correlation
   - Survival analysis
   - Multivariate modeling

### Output Structure

```
results/
├── qc/                     # Quality control reports
├── alignments/             # BAM files and indices
├── transcripts/            # Transcript identification results
├── splicing/               # Splicing analysis results
├── fusion/                 # Fusion detection results
├── ase/                    # Allele-specific expression results
├── clinical/               # Clinical integration results
└── reports/                # Summary HTML reports
```

## 📝 Configuration

The pipeline is configured through a YAML file (`config/config.yaml`):

# Reference files
reference:
  genome: "data/reference/genome.fa"
  annotation: "data/reference/annotation.gtf"
  variants: "data/reference/variants.vcf"  # Optional
```

See `config/config.example.yaml` for a complete example with all parameters.

## 🚀 Benchmarks

The pipeline has been tested on:

- ONT PromethION data (30-50M reads per sample)
- PacBio IsoSeq data (3-5M reads per sample)
- Human cancer samples (tumor/normal pairs)

Approximate runtimes on a 16-core server:
- Single sample: 4-6 hours
- Cohort of 10 samples: 24-36 hours (with parallelization)

Memory requirements:
- Peak memory usage: ~32GB for human genome analysis


## 🤝 Contributing
Contributions are welcome! Please feel free to submit a Pull Request.
1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments
- The Snakemake team for their excellent workflow management system
- Developers of FLAMES, minimap2, and other tools used in this pipeline
- The long-read sequencing community for their support and feedback
