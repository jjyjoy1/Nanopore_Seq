#!/usr/bin/env python
# train_classifier.py - Train ML model for rRNA modification classification
# Save in scripts/train_classifier.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import confusion_matrix, classification_report, roc_curve, auc
from sklearn.preprocessing import StandardScaler
import pickle
import os

# Get input and output file paths from Snakemake
data_file = snakemake.input.ml_data
output_report = snakemake.output.report
output_importance = snakemake.output.importance
output_model = snakemake.output.model

# Read the processed data
print(f"Reading processed data from {data_file}")
data = pd.read_csv(data_file)

# Separate features and target
print("Preparing data for modeling")
X = data.drop(['sample_id', 'tissue_type', 'cancer_status'], axis=1)
y = data['cancer_status']

# Get feature names for later reference
feature_names = X.columns

# Scale features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=0.25, random_state=42, stratify=y
)

print(f"Training with {X_train.shape[0]} samples, testing with {X_test.shape[0]} samples")

# Train a Random Forest classifier
print("Training Random Forest classifier")
model = RandomForestClassifier(
    n_estimators=100,
    max_depth=10,
    min_samples_split=5,
    min_samples_leaf=2,
    random_state=42,
    n_jobs=-1
)

# Perform cross-validation
cv_scores = cross_val_score(model, X_train, y_train, cv=5, scoring='accuracy')
print(f"Cross-validation accuracy: {np.mean(cv_scores):.3f} ± {np.std(cv_scores):.3f}")

# Train the final model on all training data
model.fit(X_train, y_train)

# Make predictions on the test set
y_pred = model.predict(X_test)
y_pred_proba = model.predict_proba(X_test)[:,1]

# Generate classification report
report = classification_report(y_test, y_pred)
print(report)

# Save the classification report
with open(output_report, 'w') as f:
    f.write(f"Cross-validation accuracy: {np.mean(cv_scores):.3f} ± {np.std(cv_scores):.3f}\n\n")
    f.write(report)

# Calculate feature importance
feature_importance = pd.DataFrame({
    'feature': feature_names,
    'importance': model.feature_importances_
}).sort_values('importance', ascending=False)

# Save top modifications to report
with open(output_report, 'a') as f:
    f.write("\n\nTop 15 Most Important rRNA Modification Sites:\n")
    f.write(feature_importance.head(15).to_string())

# Plot feature importance
plt.figure(figsize=(12, 8))
sns.barplot(x='importance', y='feature', data=feature_importance.head(15))
plt.title('Top 15 Most Important rRNA Modification Sites')
plt.tight_layout()
plt.savefig(output_importance, dpi=300)
print(f"Feature importance plot saved to {output_importance}")

# Calculate ROC curve and AUC
fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
roc_auc = auc(fpr, tpr)

# Plot ROC curve
plt.figure(figsize=(8, 8))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.savefig(os.path.splitext(output_importance)[0] + '_roc.png', dpi=300)

# Save the trained model
with open(output_model, 'wb') as f:
    pickle.dump(model, f)
print(f"Trained model saved to {output_model}")

# Generate confusion matrix
cm = confusion_matrix(y_test, y_pred)
plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
plt.xlabel('Predicted labels')
plt.ylabel('True labels')
plt.title('Confusion Matrix')
plt.savefig(os.path.splitext(output_importance)[0] + '_confusion.png', dpi=300)

print("Analysis complete!")

