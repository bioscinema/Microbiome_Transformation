# HPC Result Analysis

# Load necessary library
library(tidyverse)

# Read the data
data <- read.csv("all_summarized_results_GUniFrac.csv")

# Initialize list to store models
models <- list()
coefficients_list <- list()

# Loop through each transformation and fit the linear model
transformations <- unique(data$transformation)

for (trans in transformations) {
  subset_data <- data %>%
    filter(transformation == trans)
  model <- lm(mean_power ~ diff_otu_direct + diff_otu_mode + depth_mu + depth_theta + covariate_eff_sd + confounder_eff_sd + depth_conf_factor, data = subset_data)
  models[[trans]] <- model
  
  # Extract coefficients
  coefs <- summary(model)$coefficients
  coefs <- as.data.frame(coefs)
  coefs$variable <- rownames(coefs)
  coefs <- coefs %>% 
    mutate(transformation = trans) %>% 
    select(transformation, variable, Estimate)
  coefficients_list[[trans]] <- coefs
}

# Save the model summaries to a text file
output_file <- "model_summaries.txt"

# Redirect output to the file
sink(output_file)

# Print the summary of each model
for (transformation in names(models)) {
  cat("\nTransformation:", transformation, "\n")
  print(summary(models[[transformation]]))
}
sink()

cat("Model summaries have been saved to", output_file, "\n")

# Combine all coefficients into one data frame
coefficients_df <- bind_rows(coefficients_list)

# Reshape the data frame to have one row per transformation and variable names as columns
coefficients_wide <- coefficients_df %>%
  pivot_wider(names_from = variable, values_from = Estimate)

# Format all numeric values to 7 decimal places
coefficients_wide <- coefficients_wide %>%
  mutate(across(where(is.numeric), ~sprintf("%.7f", .)))

# Save the table to a CSV file
output_file_csv <- "coefficients_table.csv"
write.csv(coefficients_wide, file = output_file_csv, row.names = FALSE)

cat("Coefficients table has been saved to", output_file_csv, "\n")





models <- list()
coefficients_list <- list()

# Loop through each transformation and fit the linear model
transformations <- unique(data$transformation)

for (trans in transformations) {
  subset_data <- data %>%
    filter(transformation == trans)
  model <- lm(mean_fdr ~ diff_otu_direct + diff_otu_mode + depth_mu + depth_theta + covariate_eff_sd + confounder_eff_sd + depth_conf_factor, data = subset_data)
  models[[trans]] <- model
  
  # Extract coefficients
  coefs <- summary(model)$coefficients
  coefs <- as.data.frame(coefs)
  coefs$variable <- rownames(coefs)
  coefs <- coefs %>% 
    mutate(transformation = trans) %>% 
    dplyr::select(transformation, variable, Estimate)
  coefficients_list[[trans]] <- coefs
}

# Save the model summaries to a text file
output_file <- "model_summaries_fdr.txt"

# Redirect output to the file
sink(output_file)

# Print the summary of each model
for (transformation in names(models)) {
  cat("\nTransformation:", transformation, "\n")
  print(summary(models[[transformation]]))
}
sink()

cat("Model summaries have been saved to", output_file, "\n")

# Combine all coefficients into one data frame
coefficients_df <- bind_rows(coefficients_list)

# Reshape the data frame to have one row per transformation and variable names as columns
coefficients_wide <- coefficients_df %>%
  pivot_wider(names_from = variable, values_from = Estimate)

# Format all numeric values to 7 decimal places
coefficients_wide <- coefficients_wide %>%
  mutate(across(where(is.numeric), ~sprintf("%.7f", .)))

# Save the table to a CSV file
output_file_csv <- "coefficients_table_fdr.csv"
write.csv(coefficients_wide, file = output_file_csv, row.names = FALSE)

cat("Coefficients table has been saved to", output_file_csv, "\n")