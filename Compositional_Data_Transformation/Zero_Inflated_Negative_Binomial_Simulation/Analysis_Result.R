library(tidyverse)
models <- list()
summary_non_significant_as_reference = read.csv("summary_non_significant_as_reference.csv")
transformations <- unique(summary_non_significant_as_reference$transformation)
models <- list()
for (trans in transformations) {
  subset_data <- summary_non_significant_as_reference %>%
    filter(transformation == trans)
  model <- lm(mean_power ~ alpha + beta0 + beta + q, data = subset_data)
  models[[trans]] <- model
}
output_file <- "model_summaries_non_significant_as_reference.txt"

sink(output_file)

for (transformation in names(models)) {
  cat("\nTransformation:", transformation, "\n")
  print(summary(models[[transformation]]))
}
sink()

cat("Model summaries have been saved to", output_file, "\n")

model_summaries <- data.frame(Transformation = character(), 
                              Intercept = numeric(), 
                              alpha = numeric(), 
                              beta0 = numeric(), 
                              beta = numeric(), 
                              q = numeric(), 
                              stringsAsFactors = FALSE)

# Loop through each transformation and fit the linear regression model
for (trans in transformations) {
  subset_data <- summary_non_significant_as_reference %>%
    filter(transformation == trans)
  model <- lm(mean_power ~ alpha + beta0 + beta + q, data = subset_data)
  summary_model <- summary(model)
  
  # Extract coefficients
  intercept <- round(summary_model$coefficients[1, "Estimate"], 10)
  alpha <- round(summary_model$coefficients[2, "Estimate"], 10)
  beta0 <- round(summary_model$coefficients[3, "Estimate"], 10)
  beta <- round(summary_model$coefficients[4, "Estimate"], 10)
  q <- round(summary_model$coefficients[5, "Estimate"], 10)
  
  # Add the model summary to the data frame
  model_summaries <- rbind(model_summaries, data.frame(Transformation = trans, 
                                                       Intercept = intercept, 
                                                       alpha = alpha, 
                                                       beta0 = beta0, 
                                                       beta = beta, 
                                                       q = q))
}
write.csv(model_summaries, "model_coefficient.csv", row.names = FALSE)

cat("Model summaries have been saved to model_coefficient.csv\n")




library(tidyverse)
models <- list()
summary_non_significant_as_reference = read.csv("summary_non_significant_as_reference.csv")
transformations <- unique(summary_non_significant_as_reference$transformation)
models <- list()
for (trans in transformations) {
  subset_data <- summary_non_significant_as_reference %>%
    filter(transformation == trans)
  model <- lm(mean_fdr ~ alpha + beta0 + beta + q, data = subset_data)
  models[[trans]] <- model
}
output_file <- "model_summaries_non_significant_as_reference_fdr.txt"

sink(output_file)

for (transformation in names(models)) {
  cat("\nTransformation:", transformation, "\n")
  print(summary(models[[transformation]]))
}
sink()

cat("Model summaries have been saved to", output_file, "\n")

model_summaries <- data.frame(Transformation = character(), 
                              Intercept = numeric(), 
                              alpha = numeric(), 
                              beta0 = numeric(), 
                              beta = numeric(), 
                              q = numeric(), 
                              stringsAsFactors = FALSE)

# Loop through each transformation and fit the linear regression model
for (trans in transformations) {
  subset_data <- summary_non_significant_as_reference %>%
    filter(transformation == trans)
  model <- lm(mean_fdr ~ alpha + beta0 + beta + q, data = subset_data)
  summary_model <- summary(model)
  
  # Extract coefficients
  intercept <- round(summary_model$coefficients[1, "Estimate"], 10)
  alpha <- round(summary_model$coefficients[2, "Estimate"], 10)
  beta0 <- round(summary_model$coefficients[3, "Estimate"], 10)
  beta <- round(summary_model$coefficients[4, "Estimate"], 10)
  q <- round(summary_model$coefficients[5, "Estimate"], 10)
  
  # Add the model summary to the data frame
  model_summaries <- rbind(model_summaries, data.frame(Transformation = trans, 
                                                       Intercept = intercept, 
                                                       alpha = alpha, 
                                                       beta0 = beta0, 
                                                       beta = beta, 
                                                       q = q))
}
write.csv(model_summaries, "model_coefficient_fdr.csv", row.names = FALSE)

cat("Model summaries have been saved to model_coefficient_fdr.csv\n")