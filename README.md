Here’s your document converted into a clear, structured **Markdown (MD)** format:

---

# **Genomics Workflows Documentation**

## **1. Hybrid Nanopore-Illumina Bacterial Genome Assembly**

### **Pipeline Overview**
A Snakemake pipeline for hybrid assembly using long (Nanopore) and short (Illumina) reads.

#### **Key Steps**
1. **Quality Control**  
   - **Nanopore**: `NanoPlot` (QC) → `Filtlong` (filtering)  
   - **Illumina**: `fastp` (adapter trimming + QC)  

2. **Assembly & Polishing**  
   - **Initial Assembly**: `Flye` (long-read assembly)  
   - **Polishing**:  
     - Illumina: `Pilon` (error correction)  
     - Nanopore: `Racon` + `Medaka` (structural refinement)  

3. **Quality Assessment**  
   - `QUAST` (assembly metrics)  
   - `BUSCO` (completeness)  

#### **Software Requirements**
```plaintext
Snakemake (v7.0+), Python (v3.8+), NanoPlot, Filtlong, fastp, Flye,  
BWA, Samtools, Pilon, Minimap2, Racon, Medaka, QUAST, BUSCO
```

#### **Configuration (`config.yaml`)**
```yaml
samples:
  sample1:
    nanopore: "path/to/nanopore.fastq"
    illumina: ["path/to/illumina_R1.fastq", "path/to/illumina_R2.fastq"]
genome_size: "5m"  # Estimated genome size
busco_lineage: "bacteria_odb10"
threads: 16
```

---

## **2. Direct RNA Sequencing for Cancer Biomarkers**

### **Workflow Highlights**
- **Technology**: Nanopore direct RNA-seq (native RNA, no cDNA conversion).  
- **Applications**:  
  - Detect rRNA modifications (e.g., m6A, m5C) as cancer biomarkers.  
  - Tissue-specific "fingerprinting" via modification patterns.  

#### **Analysis Steps**
1. **Basecalling**: `Guppy` (modification-aware).  
2. **Mod Detection**: `Tombo` (signal analysis) → `Nanopolish` (event alignment).  
3. **Machine Learning**:  
   - Train classifiers to distinguish cancer vs. normal tissues.  

---

## **3. Long-Read RNA-Seq in Cancer Genomics**

### **Four Core Modules**
| **Module**               | **Tools**            | **Output**                          |
|--------------------------|----------------------|-------------------------------------|
| **Transcript Isoforms**  | FLAMES, StringTie2   | Full-length transcript quantification |
| **Splicing Analysis**    | SUPPA2, rMATS        | Tumor-specific splice junctions     |
| **Fusion Detection**     | Arriba, STAR-Fusion  | Gene fusions (e.g., EGFR-SEPT14)    |
| **Allele-Specific Exp.** | ASEP, Nanopore-Phasing | Allelic imbalance (e.g., TP53)     |

### **1. Setup**
```bash
conda env create -f envs/nanopore.yaml
snakemake --cores 32 --use-conda
```

### **2. Outputs**
- **Reports**: `MultiQC` (aggregates QC stats).  
- **Visualization**: `IGV` (validate fusions/isoforms).  

---

## **Key Advantages**
- **Reproducible**: Version-controlled via Conda/Snakemake.  
- **Scalable**: Parallel execution for multi-sample projects.  
- **Comprehensive**: QC at every step (raw reads → clinical insights).  

---

### **Need Further Customization?**
- For **publications**: I can provide LaTeX-ready methods text.  

Let me know if you'd like to refine any section!
