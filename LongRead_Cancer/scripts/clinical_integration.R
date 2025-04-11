#!/usr/bin/env Rscript
# Integrate RNA-seq results with clinical data for survival analysis

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)

# Snakemake integration
if (exists("snakemake")) {
  # Get input and output file paths from Snakemake
  isoform_files <- snakemake@input[["isoforms"]]
  junction_files <- snakemake@input[["junctions"]]
  fusion_files <- snakemake@input[["fusions"]]
  ase_files <- snakemake@input[["ase"]]
  clinical_file <- snakemake@input[["clinical"]]
  
  output_integrated <- snakemake@output[["integrated"]]
  output_survival <- snakemake@output[["survival"]]
  plots_dir <- snakemake@output[["plots_dir"]]
  
  outcome_col <- snakemake@params[["outcome_col"]]
  event_col <- snakemake@params[["event_col"]]
  sample_id_col <- snakemake@params[["sample_id_col"]]
} else {
  # Default values for testing
  isoform_files <- list.files("results/transcripts", pattern="*_isoforms.csv", full.names=TRUE)
  junction_files <- list.files("results/splicing", pattern="*_junctions.csv", full.names=TRUE)
  fusion_files <- list.files("results/fusion", pattern="*_fusions.csv", full.names=TRUE)
  ase_files <- list.files("results/ase", pattern="*_ase_results.csv", full.names=TRUE)
  clinical_file <- "data/clinical/clinical_data.csv"
  
  output_integrated <- "results/clinical/integrated_data.csv"
  output_survival <- "results/clinical/survival_analysis.csv"
  plots_dir <- "results/clinical/plots"
  
  outcome_col <- "survival_months"
  event_col <- "event"
  sample_id_col <- "sample_id"
}

# Create plots directory
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Load clinical data
clinical_data <- read.csv(clinical_file)

# Process isoform data
process_isoforms <- function(files) {
  isoform_data_list <- list()
  
  for (file in files) {
    # Extract sample ID from filename
    sample_id <- gsub("_isoforms\\.csv$", "", basename(file))
    
    # Read data
    isoform_data <- read.csv(file)
    
    if (nrow(isoform_data) > 0) {
      # Add sample ID
      isoform_data$sample_id <- sample_id
      
      # Keep only the columns we need
      isoform_data <- isoform_data %>%
        select(sample_id, gene_id, transcript_id, abundance) %>%
        # Calculate relative isoform usage within gene
        group_by(sample_id, gene_id) %>%
        mutate(rel_abundance = abundance / sum(abundance)) %>%
        ungroup()
      
      isoform_data_list[[length(isoform_data_list) + 1]] <- isoform_data
    }
  }
  
  # Combine all samples
  if (length(isoform_data_list) > 0) {
    combined_isoforms <- bind_rows(isoform_data_list)
    
    # Pivot to wide format for isoform fractions
    isoform_wide <- combined_isoforms %>%
      select(sample_id, transcript_id, rel_abundance) %>%
      pivot_wider(names_from = transcript_id, 
                  values_from = rel_abundance, 
                  names_prefix = "iso_")
    
    return(isoform_wide)
  } else {
    return(NULL)
  }
}

# Process junction data
process_junctions <- function(files) {
  junction_data_list <- list()
  
  for (file in files) {
    # Extract sample ID from filename
    sample_id <- gsub("_junctions\\.csv$", "", basename(file))
    
    # Read data
    junction_data <- read.csv(file)
    
    if (nrow(junction_data) > 0) {
      # Add sample ID
      junction_data$sample_id <- sample_id
      
      # Keep only the columns we need
      junction_data <- junction_data %>%
        select(sample_id, junction_id, read_count)
      
      junction_data_list[[length(junction_data_list) + 1]] <- junction_data
    }
  }
  
  # Combine all samples
  if (length(junction_data_list) > 0) {
    combined_junctions <- bind_rows(junction_data_list)
    
    # Pivot to wide format for junction counts
    junction_wide <- combined_junctions %>%
      pivot_wider(names_from = junction_id, 
                  values_from = read_count, 
                  names_prefix = "junc_")
    
    return(junction_wide)
  } else {
    return(NULL)
  }
}

# Process fusion data
process_fusions <- function(files) {
  fusion_data_list <- list()
  
  for (file in files) {
    # Extract sample ID from filename
    sample_id <- gsub("_fusions\\.csv$", "", basename(file))
    
    # Read data
    fusion_data <- read.csv(file)
    
    if (nrow(fusion_data) > 0) {
      # Add sample ID
      fusion_data$sample_id <- sample_id
      
      # Keep only the columns we need
      fusion_data <- fusion_data %>%
        select(sample_id, fusion_name, supporting_reads)
      
      fusion_data_list[[length(fusion_data_list) + 1]] <- fusion_data
    }
  }
  
  # Combine all samples
  if (length(fusion_data_list) > 0) {
    combined_fusions <- bind_rows(fusion_data_list)
    
    # Pivot to wide format for fusion counts
    fusion_wide <- combined_fusions %>%
      pivot_wider(names_from = fusion_name, 
                  values_from = supporting_reads, 
                  names_prefix = "fusion_")
    
    return(fusion_wide)
  } else {
    return(NULL)
  }
}

# Process ASE data
process_ase <- function(files) {
  ase_data_list <- list()
  
  for (file in files) {
    # Extract sample ID from filename
    sample_id <- gsub("_ase_results\\.csv$", "", basename(file))
    
    # Read data
    ase_data <- read.csv(file)
    
    if (nrow(ase_data) > 0) {
      # Add sample ID
      ase_data$sample_id <- sample_id
      
      # Keep only the columns we need
      ase_data <- ase_data %>%
        select(sample_id, gene_id, ase_score)
      
      ase_data_list[[length(ase_data_list) + 1]] <- ase_data
    }
  }
  
  # Combine all samples
  if (length(ase_data_list) > 0) {
    combined_ase <- bind_rows(ase_data_list)
    
    # Pivot to wide format for ASE scores
    ase_wide <- combined_ase %>%
      pivot_wider(names_from = gene_id, 
                  values_from = ase_score, 
                  names_prefix = "ase_")
    
    return(ase_wide)
  } else {
    return(NULL)
  }
}

# Process data
isoform_wide <- process_isoforms(isoform_files)
junction_wide <- process_junctions(junction_files)
fusion_wide <- process_fusions(fusion_files)
ase_wide <- process_ase(ase_files)

# Integrate with clinical data
integrated_data <- clinical_data

if (!is.null(isoform_wide)) {
  integrated_data <- integrated_data %>%
    left_join(isoform_wide, by = sample_id_col)
}

if (!is.null(junction_wide)) {
  integrated_data <- integrated_data %>%
    left_join(junction_wide, by = sample_id_col)
}

if (!is.null(fusion_wide)) {
  integrated_data <- integrated_data %>%
    left_join(fusion_wide, by = sample_id_col)
}

if (!is.null(ase_wide)) {
  integrated_data <- integrated_data %>%
    left_join(ase_wide, by = sample_id_col)
}

# Save integrated data
write.csv(integrated_data, output_integrated, row.names = FALSE)

# Perform survival analysis
survival_results <- data.frame()

# Check if clinical data has survival information
if (outcome_col %in% colnames(clinical_data) && event_col %in% colnames(clinical_data)) {
  # Create survival object
  surv_obj <- Surv(time = integrated_data[[outcome_col]], event = integrated_data[[event_col]])
  
  # Analyze associations with isoforms
  if (!is.null(isoform_wide)) {
    iso_cols <- colnames(isoform_wide)[grepl("^iso_", colnames(isoform_wide))]
    
    for (col in iso_cols) {
      # Skip if too many missing values
      if (sum(is.na(integrated_data[[col]])) > nrow(integrated_data) * 0.3) {
        next
      }
      
      # Try Cox regression
      tryCatch({
        cox_model <- coxph(surv_obj ~ integrated_data[[col]])
        
        # Extract statistics
        hr <- exp(coef(cox_model))
        ci <- exp(confint(cox_model))
        p_val <- summary(cox_model)$coefficients[5]
        
        # Add to results
        result_row <- data.frame(
          feature = col,
          hazard_ratio = hr,
          ci_lower = ci[1],
          ci_upper = ci[2],
          p_value = p_val
        )
        
        survival_results <- rbind(survival_results, result_row)
      }, error = function(e) {
        # Skip if error
      })
    }
  }
  
  # Analyze associations with junctions
  if (!is.null(junction_wide)) {
    junction_cols <- colnames(junction_wide)[grepl("^junc_", colnames(junction_wide))]
    
    for (col in junction_cols) {
      # Skip if too many missing values
      if (sum(is.na(integrated_data[[col]])) > nrow(integrated_data) * 0.3) {
        next
      }
      
      # Try Cox regression
      tryCatch({
        cox_model <- coxph(surv_obj ~ integrated_data[[col]])
        
        # Extract statistics
        hr <- exp(coef(cox_model))
        ci <- exp(confint(cox_model))
        p_val <- summary(cox_model)$coefficients[5]
        
        # Add to results
        result_row <- data.frame(
          feature = col,
          hazard_ratio = hr,
          ci_lower = ci[1],
          ci_upper = ci[2],
          p_value = p_val
        )
        
        survival_results <- rbind(survival_results, result_row)
      }, error = function(e) {
        # Skip if error
      })
    }
  }
  
  # Analyze associations with fusions
  if (!is.null(fusion_wide)) {
    fusion_cols <- colnames(fusion_wide)[grepl("^fusion_", colnames(fusion_wide))]
    
    for (col in fusion_cols) {
      # Skip if too many missing values
      if (sum(is.na(integrated_data[[col]])) > nrow(integrated_data) * 0.3) {
        next
      }
      
      # Try Cox regression
      tryCatch({
        cox_model <- coxph(surv_obj ~ integrated_data[[col]])
        
        # Extract statistics
        hr <- exp(coef(cox_model))
        ci <- exp(confint(cox_model))
        p_val <- summary(cox_model)$coefficients[5]
        
        # Add to results
        result_row <- data.frame(
          feature = col,
          hazard_ratio = hr,
          ci_lower = ci[1],
          ci_upper = ci[2],
          p_value = p_val
        )
        
        survival_results <- rbind(survival_results, result_row)
      }, error = function(e) {
        # Skip if error
      })
    }
  }
  
  # Analyze associations with ASE
  if (!is.null(ase_wide)) {
    ase_cols <- colnames(ase_wide)[grepl("^ase_", colnames(ase_wide))]
    
    for (col in ase_cols) {
      # Skip if too many missing values
      if (sum(is.na(integrated_data[[col]])) > nrow(integrated_data) * 0.3) {
        next
      }
      
      # Try Cox regression
      tryCatch({
        cox_model <- coxph(surv_obj ~ integrated_data[[col]])
        
        # Extract statistics
        hr <- exp(coef(cox_model))
        ci <- exp(confint(cox_model))
        p_val <- summary(cox_model)$coefficients[5]
        
        # Add to results
        result_row <- data.frame(
          feature = col,
          hazard_ratio = hr,
          ci_lower = ci[1],
          ci_upper = ci[2],
          p_value = p_val
        )
        
        survival_results <- rbind(survival_results, result_row)
      }, error = function(e) {
        # Skip if error
      })
    }
  }
  
  # Create Kaplan-Meier plots for top significant features
  if (nrow(survival_results) > 0) {
    # Adjust p-values
    survival_results$padj <- p.adjust(survival_results$p_value, method = "BH")
    
    # Sort by adjusted p-value
    survival_results <- survival_results %>%
      arrange(padj)
    
    # Create KM plots for top 5 significant features
    top_features <- head(survival_results, 5)
    
    for (i in 1:nrow(top_features)) {
      feature <- top_features$feature[i]
      
      # Create high/low groups based on median
      threshold <- median(integrated_data[[feature]], na.rm = TRUE)
      integrated_data$group <- ifelse(integrated_data[[feature]] > threshold, "High", "Low")
      
      # Create KM plot
      km_fit <- survfit(surv_obj ~ group, data = integrated_data)
      
      # Plot
      km_plot <- ggsurvplot(
        km_fit,
        data = integrated_data,
        pval = TRUE,
        risk.table = TRUE,
        title = paste0("Survival by ", feature),
        xlab = "Time (months)",
        legend.labs = c("Low", "High"),
        legend.title = feature,
        palette = c("#2E9FDF", "#E7B800")
      )
      
      # Save plot
      pdf(file.path(plots_dir, paste0("survival_", gsub("[::]", "_", feature), ".pdf")), 
          width = 10, height = 8)
      print(km_plot)
      dev.off()
    }
  }
  
  # Create multivariate models for key features
  if (nrow(survival_results) > 0) {
    # Select top features from each category for multivariate model
    top_features_by_category <- data.frame()
    
    # Top isoform
    if (!is.null(isoform_wide)) {
      iso_results <- survival_results %>%
        filter(grepl("^iso_", feature)) %>%
        arrange(padj) %>%
        head(1)
      top_features_by_category <- bind_rows(top_features_by_category, iso_results)
    }
    
    # Top junction
    if (!is.null(junction_wide)) {
      junc_results <- survival_results %>%
        filter(grepl("^junc_", feature)) %>%
        arrange(padj) %>%
        head(1)
      top_features_by_category <- bind_rows(top_features_by_category, junc_results)
    }
    
    # Top fusion
    if (!is.null(fusion_wide)) {
      fusion_results <- survival_results %>%
        filter(grepl("^fusion_", feature)) %>%
        arrange(padj) %>%
        head(1)
      top_features_by_category <- bind_rows(top_features_by_category, fusion_results)
    }
    
    # Top ASE
    if (!is.null(ase_wide)) {
      ase_results <- survival_results %>%
        filter(grepl("^ase_", feature)) %>%
        arrange(padj) %>%
        head(1)
      top_features_by_category <- bind_rows(top_features_by_category, ase_results)
    }
    
    # Perform multivariate analysis if we have at least 2 features
    if (nrow(top_features_by_category) >= 2) {
      # Create formula for multivariate model
      features <- top_features_by_category$feature
      formula_str <- paste0("surv_obj ~ ", paste(paste0("integrated_data$", features), collapse = " + "))
      formula_obj <- as.formula(formula_str)
      
      # Create multivariate Cox model
      tryCatch({
        multi_cox <- coxph(formula_obj)
        
        # Save summary to file
        sink(file.path(plots_dir, "multivariate_model.txt"))
        print(summary(multi_cox))
        sink()
        
        # Create forest plot for multivariate model
        forest_data <- as.data.frame(summary(multi_cox)$coefficients)
        forest_data$feature <- rownames(forest_data)
        forest_data$hazard_ratio <- exp(forest_data$coef)
        forest_data$ci_lower <- exp(coef(multi_cox) - 1.96 * sqrt(diag(vcov(multi_cox))))
        forest_data$ci_upper <- exp(coef(multi_cox) + 1.96 * sqrt(diag(vcov(multi_cox))))
        
        # Create forest plot
        p <- ggplot(forest_data, aes(x = hazard_ratio, y = feature)) +
          geom_point(size = 3) +
          geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
          geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
          scale_x_log10() +
          labs(x = "Hazard Ratio (log scale)", y = "", title = "Multivariate Survival Analysis") +
          theme_bw(base_size = 12)
        
        # Save plot
        ggsave(file.path(plots_dir, "multivariate_forest_plot.pdf"), p, width = 10, height = 6)
      }, error = function(e) {
        # Print error message if multivariate model fails
        message("Failed to create multivariate model: ", e$message)
      })
    }
  }
}

# Save survival results
write.csv(survival_results, output_survival, row.names = FALSE)

# Create summary plot for survival associations
if (nrow(survival_results) > 0) {
  # Filter for significant results
  sig_results <- survival_results %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(log(hazard_ratio))))
  
  if (nrow(sig_results) > 0) {
    # Top 20 features
    top_sig <- head(sig_results, 20)
    
    # Create forest plot
    p <- ggplot(top_sig, aes(x = hazard_ratio, y = reorder(feature, hazard_ratio))) +
      geom_point(size = 3) +
      geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      scale_x_log10() +
      labs(x = "Hazard Ratio (log scale)", y = "", title = "Survival Associations") +
      theme_bw(base_size = 12) +
      theme(axis.text.y = element_text(size = 10))
    
    # Save plot
    ggsave(file.path(plots_dir, "survival_forest_plot.pdf"), p, width = 10, height = 8)
    
    # Create table with feature categories
    feature_categories <- data.frame(
      feature = top_sig$feature,
      category = ifelse(grepl("^iso_", top_sig$feature), "Isoform",
                ifelse(grepl("^junc_", top_sig$feature), "Splice Junction",
                ifelse(grepl("^fusion_", top_sig$feature), "Fusion Transcript",
                ifelse(grepl("^ase_", top_sig$feature), "Allele-Specific Expression", "Other")))),
      hazard_ratio = top_sig$hazard_ratio,
      p_value = top_sig$p_value,
      padj = top_sig$padj
    )
    
    # Save table
    write.csv(feature_categories, file.path(plots_dir, "survival_significant_features.csv"), row.names = FALSE)
    
    # Create barplot by category
    category_counts <- feature_categories %>%
      count(category) %>%
      arrange(desc(n))
    
    p_cat <- ggplot(category_counts, aes(x = reorder(category, n), y = n, fill = category)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(x = "", y = "Number of Significant Features", title = "Significant Features by Category") +
      theme_minimal()
    
    # Save plot
    ggsave(file.path(plots_dir, "survival_feature_categories.pdf"), p_cat, width = 8, height = 6)
  }
}

# Print completion message
cat("Clinical integration and survival analysis completed.\n")

