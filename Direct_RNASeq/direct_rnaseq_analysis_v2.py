# Step 1: Start with raw FAST5 files from the sequencer
# These contain the electrical signals needed for modification detection

# Step 2: Base-calling (if not already done by the sequencer)
# Using Guppy or other ONT base-callers
os.system("guppy_basecaller -i fast5_files/ -s basecalled/ --config rna_r9.4.1_70bps_hac.cfg --device cuda:0")

# Step 3: Align reads to reference transcriptome
# Using minimap2, which is optimized for nanopore data
os.system("minimap2 -ax map-ont reference_transcriptome.fa basecalled/pass/*.fastq > alignments.sam")
os.system("samtools sort -o alignments.sorted.bam alignments.sam")
os.system("samtools index alignments.sorted.bam")

# Step 4: Resquiggle - crucial step that realigns raw signal to reference sequence
os.system("tombo resquiggle fast5_files/ reference_transcriptome.fa --processes 8 --fit-scale --include-event-stdev --corrected-group RescoreFast5")

# Step 5: Detect modifications based on signal deviations
os.system("tombo detect_modifications alternative_model --fast5-basedirs fast5_files/ " + 
          "--statistics-file-basename mod_stats " + 
          "--alternate-bases all")

# Step 6: Extract modification data to a format usable for machine learning
# This would generate a per-position modification probability table
os.system("tombo text_output dampened_fraction --statistics-filename mod_stats.rna.tombo.stats " + 
          "--output-basename rRNA_mods " + 
          "--region-type transcriptome")

# Now we'd have a file like 'rRNA_mods.dampened_fraction.csv' with modification probabilities
# Format: ref_id, pos, fraction_modified

# Convert this to a sample-by-modification matrix for machine learning
# This would require custom code to reshape the data:

import pandas as pd

# Read modification data
mod_data = pd.read_csv('rRNA_mods.dampened_fraction.csv', sep='\t')

# Create a pivot table: rows=samples, columns=modification sites
mod_matrix = mod_data.pivot_table(
    index='sample_id',  # You'd need to add this column during processing
    columns=['ref_id', 'pos'],
    values='fraction_modified',
    fill_value=0
)

# Add metadata
metadata = pd.read_csv('sample_metadata.csv')  # Contains tissue_type, cancer_status
full_dataset = pd.merge(metadata, mod_matrix, on='sample_id')

# Save for machine learning
full_dataset.to_csv('rRNA_modifications_for_ML.csv')
