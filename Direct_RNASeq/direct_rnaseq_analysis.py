# Installation of key packages
# pip install ont-tombo numpy pandas matplotlib scikit-learn seaborn

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report

# Step 1: Preprocess the raw nanopore data using Tombo
# Assuming data is in 'fast5_files' directory
os.system("tombo resquiggle fast5_files reference_transcriptome.fa --processes 8 --fit-scale --include-event-stdev")

# Step 2: Detect modified positions using Tombo
os.system("tombo detect_modifications alternative_model --fast5-basedirs fast5_files " + 
          "--statistics-file-basename mod_stats " + 
          "--alternate-bases all")

# Step 3: Extract features from the modification detection results
os.system("tombo text_output browser_files --fast5-basedirs fast5_files " + 
          "--statistics-filename mod_stats.rna.5mC.tombo.stats " + 
          "--browser-file-basename browser_files")

# Step 4: Load and prepare the data for machine learning
# This assumes you've extracted features into a CSV with columns for each modification site
modifications_df = pd.read_csv('rRNA_modifications.csv')

# Assume columns are: sample_id, tissue_type, cancer_status, mod_site_1, mod_site_2, etc.
X = modifications_df.iloc[:, 3:]  # All modification sites
y = modifications_df['cancer_status']  # Target variable

# Step 5: Train a classifier
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42)

model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Step 6: Evaluate the model
predictions = model.predict(X_test)
print(classification_report(y_test, predictions))

# Step 7: Visualize important modification sites
feature_importance = pd.DataFrame({
    'feature': X.columns,
    'importance': model.feature_importances_
}).sort_values('importance', ascending=False)

plt.figure(figsize=(10, 8))
sns.barplot(x='importance', y='feature', data=feature_importance.head(15))
plt.title('Top 15 Most Important rRNA Modification Sites')
plt.tight_layout()
plt.savefig('important_modifications.png')

