library(tidyverse)
models <- list()
summary_significant_as_reference = read.csv("summary_significant_as_reference.csv")

transformations <- unique(summary_significant_as_reference$transformation)

for (trans in transformations) {
  subset_data <- summary_significant_as_reference %>%
    filter(transformation == trans)
  model <- lm(mean_power ~ alpha + beta0 + beta + q, data = subset_data)
  models[[trans]] <- model
}
output_file <- "model_summaries_significant_as_reference.txt"

sink(output_file)

for (transformation in names(models)) {
  cat("\nTransformation:", transformation, "\n")
  print(summary(models[[transformation]]))
}
sink()

cat("Model summaries have been saved to", output_file, "\n")

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

