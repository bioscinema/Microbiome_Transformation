library(phyloseq) 
library(DESeq2)
library(dplyr)
library(tibble) 
library(tidyverse)
library(reticulate)
load("IBD_16s_data_V3.RData")
replace_zeros <- function(df, constant = 1) {
  df[df == 0] <- constant
  return(df)
}
filter_rows <- function(data) {
  filtered_data <- data[apply(data, 1, function(row) {
    sum(row != 0) >= 20
  }), ]
  return(filtered_data)
}
otu_data <- data.frame(replace_zeros(filter_rows(otu_table(phy1))))
sample_metadata <- sample_data(phy1)
sample_metadata[["diagnosis"]] = ifelse(sample_metadata[["diagnosis"]] == "healthy_control", "A", "B")
rownames(sample_metadata) = paste0("X", rownames(sample_metadata))

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(otu_data),
  colData = sample_metadata,
  design = ~ diagnosis
)
#dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds)
res <- results(dds)
significant_results <- res %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  filter(padj < 0.05)
deseq2_significant_replace1 <- as.vector(significant_results["OTU"])$OTU
deseq2_significant_replace1

otu_data <- data.frame((filter_rows(otu_table(phy1))))
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(otu_data),
  colData = sample_metadata,
  design = ~ diagnosis
)

dds <- estimateSizeFactors(dds, type="poscount")
dds <- DESeq(dds)
res <- results(dds)

# Extract significant results
significant_results <- res %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  filter(padj < 0.05)


# View the significant results
deseq2_significant_poscount <- as.vector(significant_results["OTU"])$OTU
deseq2_significant_poscount




alr_transformation <- function(data, group_factor = 1, ref_component = ncol(data)-1) {
  data = replace_zeros(data)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  factors = as.vector(unique(data[[group_factor]]))
  group1 = data[data[group_factor] == factors[1],-group_factor]
  group2 = data[data[group_factor] == factors[2],-group_factor]
  group1_alr <- log(group1[, -ref_component] / group1[, ref_component])
  group2_alr <- log(group2[, -ref_component] / group2[, ref_component])
  colnames(group1_alr) <- colnames(data[,-c(group_factor,(ref_component+1))])
  colnames(group2_alr) <- colnames(data[,-c(group_factor,(ref_component+1))])
  group1_alr <- cbind(Group = factors[1], group1_alr)
  group2_alr <- cbind(Group = factors[2], group2_alr)
  combined <- rbind(group1_alr, group2_alr)
  return(combined)
}


logit <- function(x){
  log(x/(1-x))
}

logit_alr_transformation <- function(data, group_factor = 1, ref_component = ncol(data)-1) {
  data = replace_zeros(data)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  factors = as.vector(unique(data[[group_factor]]))
  group1 = data[data[group_factor] == factors[1],-group_factor]
  group2 = data[data[group_factor] == factors[2],-group_factor]
  group1_alr <- logit(group1[, -ref_component]) - logit(group1[, ref_component])
  group2_alr <- logit(group2[, -ref_component]) - logit(group2[, ref_component])
  group1_alr <- cbind(Group = factors[1], group1_alr)
  group2_alr <- cbind(Group = factors[2], group2_alr)
  combined <- rbind(group1_alr, group2_alr)
  return(combined)
}

arcsine <- function(x) {
  return(asin(sqrt(x)))
}
arcsine_alr_transformation <- function(data, group_factor = 1) {
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  data_without_group = as.matrix(data[,-group_factor])
  p = dim(data_without_group)[2]
  T = rbind(diag(p-1),rep(-1,p-1))
  t_data = arcsine(data_without_group) %*% T
  transformed_data = data.frame(data[,group_factor], t_data)
  colnames(transformed_data)[1] <- "Group"
  colnames(transformed_data)[-1] <- colnames(data_without_group[,-ncol(data_without_group)])
  return(transformed_data)
}

clr_transformation <- function(data, group_factor = 1) {
  data = replace_zeros(data)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  factors = as.vector(unique(data[[group_factor]]))
  group1 = data[data[group_factor] == factors[1],-group_factor]
  group2 = data[data[group_factor] == factors[2],-group_factor]
  gm_group1 <- exp(rowMeans(log(group1)))
  gm_group2 <- exp(rowMeans(log(group2)))
  clr_group1 <- log(group1 / gm_group1)
  clr_group2 <- log(group2 / gm_group2)
  clr_group1 <- cbind(Group = factors[1], clr_group1)
  clr_group2 <- cbind(Group = factors[2], clr_group2)
  combined <- rbind(clr_group1, clr_group2)
  return(combined)
}


logit_clr_transformation <- function(data, group_factor = 1) {
  data = replace_zeros(data)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  factors = as.vector(unique(data[[group_factor]]))
  group1 = data[data[group_factor] == factors[1],-group_factor]
  group2 = data[data[group_factor] == factors[2],-group_factor]
  gm_group1 <- exp(rowMeans(logit(group1)))
  gm_group2 <- exp(rowMeans(logit(group2)))
  clr_group1 <- logit(group1)- logit(gm_group1)
  clr_group2 <- logit(group2)- logit(gm_group2)
  clr_group1 <- cbind(Group = factors[1], clr_group1)
  clr_group2 <- cbind(Group = factors[2], clr_group2)
  combined <- rbind(clr_group1, clr_group2)
  return(combined)
}

arcsine <- function(x) {
  return(asin(sqrt(x)))
}
arcsine_clr_transformation <- function(data, group_factor = 1) {
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  data_without_group = as.matrix(data[,-group_factor])
  p = dim(data_without_group)[2]
  T = diag(p)-1/p*matrix(1,p,p)
  t_data = arcsine(data_without_group) %*% T  
  transformed_data = data.frame(data[,group_factor], t_data)
  colnames(transformed_data)[1] <- "Group"
  colnames(transformed_data)[-1] <- colnames(data_without_group)
  return(transformed_data)
}

library(compositions) 
ilr_transformation <- function(data, group_factor = 1) {
  data = replace_zeros(data)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  pure_data = data[, -group_factor, drop = FALSE]
  ilr_data = ilr(pure_data)
  factor_col = data[, group_factor, drop = FALSE]
  names(factor_col) <- "Group"
  combined = data.frame(factor_col, ilr_data)
  colnames(combined)[-1] <- paste0("X", 1:ncol(ilr_data))
  return(combined)
}

source_python("Scource_Code.py")

dual_group_boxcox_transformation <- function(data, group_factor = 1){
  factors = as.vector(unique(data[[group_factor]]))
  n = nrow(data)
  group1 = data[data[group_factor] == factors[1],-group_factor]
  group2 = data[data[group_factor] == factors[2],-group_factor]
  transformed = matrix(nrow = nrow(data), ncol = ncol(data)-1)
  for(i in 1:((ncol(data)-1))){
    
    py$p1 = group1[,i]
    py$p2 = group2[,i]
    
    transformed_results <- py$dual_group_boxcox_transformation2(py$p1, py$p2)
    transformed[1:nrow(group1), i] <- transformed_results[[1]]
    transformed[(nrow(group1) + 1):n, i] <- transformed_results[[2]]
  }
  group_column = c(rep(factors[1], nrow(group1)), rep(factors[2], nrow(group2)))
  transformed_df <- data.frame(Group = group_column, transformed)
  return(data.frame(transformed_df))
}

dual_group_boxcox_alr_transformation <- function(data, group_factor = 1){
  data = replace_zeros(data)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  boxcox_transformed = dual_group_boxcox_transformation(data, group_factor = 1)
  data_without_group = as.matrix(boxcox_transformed[,-group_factor])
  p = dim(data_without_group)[2]
  T = rbind(diag(p-1),rep(-1,p-1))
  t_data = (data_without_group) %*% T
  transformed_data = data.frame(data[,group_factor], t_data)
  colnames(transformed_data)[1] <- "Group"
  colnames(transformed_data)[-1] <- colnames(data_without_group[,-ncol(data_without_group)])
  return(transformed_data)  
}

dual_group_boxcox_clr_transformation <- function(data, group_factor = 1){
  data = replace_zeros(data)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  boxcox_transformed = dual_group_boxcox_transformation(data, group_factor = 1)
  data_without_group = as.matrix(boxcox_transformed[,-group_factor])
  p = dim(data_without_group)[2]
  T = diag(p)-1/p*matrix(1,p,p)
  t_data = (data_without_group) %*% T
  transformed_data = data.frame(data[,group_factor], t_data)
  colnames(transformed_data)[1] <- "Group"
  colnames(transformed_data)[-1] <- colnames(data_without_group)
  return(transformed_data)  
}

library(MASS)
boxcoxt <- function(data) {
  if(any(data <= 0)) {
    stop("Data must contain only positive values for Box-Cox transformation.")
  }
  bc_out <- boxcox(data ~ 1, plotit = FALSE)
  lambda_optimal <- bc_out$x[which.max(bc_out$y)]
  # Apply the Box-Cox transformation
  if(lambda_optimal == 0) {
    transformed_data <- log(data)
  } else {
    transformed_data <- (data^lambda_optimal - 1) / lambda_optimal
  }
  
  return(transformed_data)
}
boxcox_alr_transformation_william <- function(data, group_factor = 1, ref_component = ncol(data)-1){
  data = replace_zeros(data)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  groups <- data[[group_factor]]
  pure_data <- data[,-group_factor]
  ratios <- pure_data[, -ref_component] / pure_data[, ref_component]
  transformed_ratios <- apply(ratios, 2, boxcoxt)
  transformed_data <- data.frame(Group = groups, transformed_ratios)
  return(transformed_data)
}

boxcox_alr_transformation_dual_group <- function(data, group_factor = 1, ref_component = ncol(data)-1){
  data = replace_zeros(data)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  groups <- data[[group_factor]]
  pure_data <- data[,-group_factor]
  ratios <- (pure_data[, -ref_component]) / pure_data[, ref_component]
  ratio_with_group = data.frame(Group = groups, ratios)
  transformed_data = dual_group_boxcox_transformation(ratio_with_group)
  return(transformed_data)
}

boxcox_transform <- function(data, group_factor = 1) {
  data = replace_zeros(data)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  n = nrow(data)
  factors = as.vector(unique(data[[group_factor]]))
  data_without_factor = data[,-group_factor]
  transformed = matrix(nrow = nrow(data), ncol = ncol(data)-1)
  for(i in 1:((ncol(data)-1))){
    col_data = data_without_factor[,i]
    transformed_data = boxcoxt(col_data)
    transformed[1:nrow(data), i] <- transformed_data
  }
  group_column = c(rep(factors[1], nrow(data[data[group_factor] == factors[1],-group_factor])), rep(factors[2], nrow(data[data[group_factor] == factors[2],-group_factor])))
  transformed_df <- data.frame(Group = group_column, transformed)
  colnames(transformed_df) = colnames(data)
  return(transformed_df)
}


power_alr_transformation <- function(data, group_factor = 1){
  boxcox_transformed = boxcox_transform(data, group_factor = 1)
  data_without_group = as.matrix(boxcox_transformed[,-group_factor])
  p = dim(data_without_group)[2]
  T = rbind(diag(p-1),rep(-1,p-1))
  t_data = (data_without_group) %*% T
  transformed_data = data.frame(data[,group_factor], t_data)
  colnames(transformed_data)[1] <- "Group"
  colnames(transformed_data)[-1] <- colnames(data_without_group[,-ncol(data_without_group)])
  return(transformed_data)  
}

power_clr_transformation <- function(data, group_factor = 1){
  boxcox_transformed = boxcox_transform(data, group_factor = 1)
  data_without_group = as.matrix(boxcox_transformed[,-group_factor])
  p = dim(data_without_group)[2]
  T = diag(p)-1/p*matrix(1,p,p)
  t_data = (data_without_group) %*% T
  transformed_data = data.frame(data[,group_factor], t_data)
  colnames(transformed_data)[1] <- "Group"
  colnames(transformed_data)[-1] <- colnames(data_without_group)
  return(transformed_data)  
}

new_logit_transformation <- function(data, group_factor = 1){
  factors = as.vector(unique(data[[group_factor]]))
  n = nrow(data)
  group1 = data[data[group_factor] == factors[1],-group_factor]
  group2 = data[data[group_factor] == factors[2],-group_factor]
  transformed = matrix(nrow = nrow(data), ncol = ncol(data)-1)
  for(i in 1:((ncol(data)-1))){
    
    py$p1 = as.numeric(group1[, i])
    py$p2 = as.numeric(group2[, i])
    
    transformed_results <- py$new_logit_transformation_s(py$p1, py$p2)
    transformed[1:nrow(group1), i] <- transformed_results[[1]]
    transformed[(nrow(group1) + 1):n, i] <- transformed_results[[2]]
  }
  group_column = c(rep(factors[1], nrow(group1)), rep(factors[2], nrow(group2)))
  transformed_df <- data.frame(Group = group_column, transformed)
  colnames(transformed_df) = colnames(data)
  return(transformed_df)
}

new_logit_alr_transformation <- function(data, group_factor = 1){
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  boxcox_transformed = new_logit_transformation(data, group_factor = 1)
  data_without_group = as.matrix(boxcox_transformed[,-group_factor])
  p = dim(data_without_group)[2]
  T = rbind(diag(p-1),rep(-1,p-1))
  t_data = (data_without_group) %*% T
  transformed_data = data.frame(data[,group_factor], t_data)
  colnames(transformed_data)[1] <- "Group"
  colnames(transformed_data)[-1] <- colnames(data_without_group[,-ncol(data_without_group)])
  return(transformed_data)  
}

new_logit_clr_transformation <- function(data, group_factor = 1){
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  boxcox_transformed = new_logit_transformation(data, group_factor = 1)
  data_without_group = as.matrix(boxcox_transformed[,-group_factor])
  p = dim(data_without_group)[2]
  T = diag(p)-1/p*matrix(1,p,p)
  t_data = (data_without_group) %*% T
  transformed_data = data.frame(data[,group_factor], t_data)
  colnames(transformed_data)[1] <- "Group"
  colnames(transformed_data)[-1] <- colnames(data_without_group)
  return(transformed_data)  
}


library(dplyr)
otu_data <- data.frame(filter_rows(otu_table(phy1)))
t_otu_data = data.frame(t(otu_data))
t_otu_data$Group <- sample_metadata$diagnosis
t_otu_data = t_otu_data %>%
  dplyr::select(Group, everything())

otu_data <- data.frame(filter_rows(otu_table(phy1)))
totu_data = t(otu_data)
totu_data = data.frame("Group" = sample_metadata[["diagnosis"]], totu_data)
generated_data <- totu_data

# Apply Wilcoxon test for each column
wilcoxon_results <- sapply(colnames(generated_data)[-1], function(col) {
  wilcox.test(generated_data[[col]] ~ generated_data$Group)$p.value
})

# Apply Mann-Whitney U test for each column
mann_whitney_results <- sapply(colnames(generated_data)[-1], function(col) {
  wilcox.test(generated_data[[col]] ~ generated_data$Group, exact = FALSE)$p.value
})

# Apply Kruskal-Wallis test for each column
kruskal_results <- sapply(colnames(generated_data)[-1], function(col) {
  kruskal.test(generated_data[[col]] ~ generated_data$Group)$p.value
})

# Apply permutation test for each column


# Apply Kolmogorov-Smirnov test for each column
ks_results <- sapply(colnames(generated_data)[-1], function(col) {
  ks.test(generated_data[[col]][generated_data$Group == "A"], generated_data[[col]][generated_data$Group == "B"])$p.value
})

# Combine results into a dataframe
test_results <- data.frame(
  Column = colnames(generated_data)[-1],
  Wilcoxon = wilcoxon_results,
  Mann_Whitney = mann_whitney_results,
  Kruskal = kruskal_results,
  Kolmogorov_Smirnov = ks_results
)

# Identify the column with the highest p-value in any test
max_p_values <- apply(test_results[-1], 1, max)
most_non_significant_column <- test_results$Column[which.max(max_p_values)]

# Print the most non-significant column and its p-values
most_non_significant_results <- test_results[test_results$Column == most_non_significant_column, ]

# Result list
results <- list(
  most_non_significant_column = most_non_significant_column,
  most_non_significant_results = most_non_significant_results,
  test_results = test_results
)

col_names <- colnames(t_otu_data)

# Identify the columns to switch
col_196 <- col_names[196]
col_last <- col_names[length(col_names)]

# Reorder the columns to switch the 196th and the last columns
col_names_switched <- col_names
col_names_switched[196] <- col_last
col_names_switched[length(col_names)] <- col_196

# Reorder the dataframe
t_otu_data <- t_otu_data[, col_names_switched]


transformation_functions <- list(
  "ALR" = alr_transformation,
  "CLR" = clr_transformation,
  "Additive Power Contrast" = power_alr_transformation,
  "Centered Power Contrast" = power_clr_transformation,
  "Additive Logit Contrast" = logit_alr_transformation,
  "Centered Logit Contrast" = logit_clr_transformation,
  "Additive Arcsine Contrast" = arcsine_alr_transformation,
  "Centered Arcsine Contrast" = arcsine_clr_transformation,
  "Additive Dual Group Logit Contrast" = new_logit_alr_transformation,
  "Centered Dual Group Logit Contrast" = new_logit_clr_transformation,
  "Additive Dual Group Power Contrast" = dual_group_boxcox_alr_transformation,
  "Centered Dual Group Power Contrast" = dual_group_boxcox_clr_transformation,
  "ILR" = ilr_transformation,
  "Boxcox in Ratio" = boxcox_alr_transformation_william,
  "Dual Group Boxcox in Ratio" = boxcox_alr_transformation_dual_group
  
)

# Function to perform Wilcoxon test and adjust p-values

# Function to perform two-sample t-test and adjust p-values
perform_t_test <- function(data) {
  # Perform t-tests for each column except the first (which is assumed to be the group variable)
  p_values <- sapply(data[-1], function(column) {
    # Check if data is essentially constant
    if (var(column[data$Group == unique(data$Group)[1]]) == 0 || var(column[data$Group == unique(data$Group)[2]]) == 0) {
      return(NA) # Return NA if data is essentially constant
    } else {
      return(t.test(column ~ data$Group)$p.value)
    }
  })
  
  # Adjust p-values using the BH method, handling NAs appropriately
  adjusted_p_values <-p.adjust(p_values, method = "BH", n = sum(!is.na(p_values)))
  
  return(adjusted_p_values)
}


# Initialize a list to store results for each transformation
results <- list()

# Loop through each transformation, apply it, and calculate the overlap and differences
for (transformation_name in names(transformation_functions)) {
  
  # Apply the transformation
  transformed_data <- transformation_functions[[transformation_name]](t_otu_data)
  
  # Perform the t-test on the transformed data
  t_test_results <- perform_t_test(transformed_data)
  
  # Determine significance based on adjusted p-value < 0.05
  t_test_significant <- colnames(transformed_data)[t_test_results < 0.05]
  
  # Calculate overlap
  overlap <- length(intersect(deseq2_significant_replace1, t_test_significant))
  
  # Calculate differences
  deseq2_only <- length(setdiff(deseq2_significant_replace1, t_test_significant))
  t_test_only <- length(setdiff(t_test_significant, deseq2_significant_replace1))
  
  # Store the results
  results[[transformation_name]] <- list(
    overlap = overlap,
    deseq2_only = deseq2_only,
    t_test_only = t_test_only
  )
}

# Convert results to a dataframe
results_df <- do.call(rbind, lapply(names(results), function(name) {
  c(Transformation = name, results[[name]])
}))

# Convert the list to a data frame
results_df <- as.data.frame(results_df, stringsAsFactors = FALSE)

# Convert numerical columns to numeric type
results_df$overlap <- as.numeric(results_df$overlap)
results_df$deseq2_only <- as.numeric(results_df$deseq2_only)
results_df$t_test_only <- as.numeric(results_df$t_test_only)

# Rename the columns for clarity
colnames(results_df) <- c("Transformation", "Overlap", "DESeq2_Only", "T_Test_Only")
results_df[,1] = unlist(results_df[,1])
# Display the results dataframe
print(results_df)

write.csv(results_df,"overlap_replace1.csv", row.names = FALSE)



# Initialize a list to store results for each transformation
results <- list()

# Loop through each transformation, apply it, and calculate the overlap and differences
for (transformation_name in names(transformation_functions)) {
  
  # Apply the transformation
  transformed_data <- transformation_functions[[transformation_name]](t_otu_data)
  
  # Perform the t-test on the transformed data
  t_test_results <- perform_t_test(transformed_data)
  
  # Determine significance based on adjusted p-value < 0.05
  t_test_significant <- colnames(transformed_data)[t_test_results < 0.05]
  
  # Calculate overlap
  overlap <- length(intersect(deseq2_significant_poscount, t_test_significant))
  
  # Calculate differences
  deseq2_only <- length(setdiff(deseq2_significant_poscount, t_test_significant))
  t_test_only <- length(setdiff(t_test_significant, deseq2_significant_poscount))
  
  # Store the results
  results[[transformation_name]] <- list(
    overlap = overlap,
    deseq2_only = deseq2_only,
    t_test_only = t_test_only
  )
}

# Convert results to a dataframe
results_df <- do.call(rbind, lapply(names(results), function(name) {
  c(Transformation = name, results[[name]])
}))

# Convert the list to a data frame
results_df <- as.data.frame(results_df, stringsAsFactors = FALSE)

# Convert numerical columns to numeric type
results_df$overlap <- as.numeric(results_df$overlap)
results_df$deseq2_only <- as.numeric(results_df$deseq2_only)
results_df$t_test_only <- as.numeric(results_df$t_test_only)

# Rename the columns for clarity
colnames(results_df) <- c("Transformation", "Overlap", "DESeq2_Only", "T_Test_Only")
results_df[,1] = unlist(results_df[,1])
# Display the results dataframe
print(results_df)

write.csv(results_df,"overlap_poscount.csv", row.names = FALSE)