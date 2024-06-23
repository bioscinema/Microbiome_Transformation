user_lib <- "/home/yxz3116/R_Package"


if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
}
chooseCRANmirror(graphics = FALSE, ind = 1)

is_installed <- function(pkg, lib) {
  return(pkg %in% rownames(installed.packages(lib.loc = lib)))
}

# Function to install and load a package if not already installed
install_and_load <- function(pkg, lib) {
  if (!is_installed(pkg, lib)) {
    install.packages(pkg, lib = lib)
  }
  library(pkg, character.only = TRUE, lib.loc = lib)
}

# List of packages to install and load
packages <- c("compositions", "reticulate", "MASS", "dplyr", "tibble")

# Install and load each package
for (pkg in packages) {
  install_and_load(pkg, user_lib)
}

generate_data <- function(n = 100, p = 500, beta0 = 5, beta = 1, alpha = 3, q = 0.3) {
  n1 <- n / 2
  n2 <- n / 2
  x <- c(rep(1, n1), rep(0, n2))
  
  # Significant columns
  mu_nonzero <- exp(beta0 + x * beta)
  size_nonzero <- 1 / alpha
  prob_nonzero <- 1 / (alpha * mu_nonzero + 1)
  xn <- replicate(p / 2, rnbinom(n, size_nonzero, prob_nonzero))
  
  # Non-significant columns
  beta_zero <- 0
  mu_zero <- exp(beta0 + x * beta_zero)
  size_zero <- 1 / alpha
  prob_zero <- 1 / (alpha * mu_zero + 1)
  xz <- replicate(p / 2, rnbinom(n, size_zero, prob_zero))
  
  X <- cbind(xn, xz)
  
  pi <- replicate(p, 1 - rbinom(n, 1, q))
  Xc <- matrix(, n, p)
  for (j in 1:p) {
    Xc[, j] <- X[, j] * pi[, j]
  }
  
  s <- rlnorm(n, meanlog = 1)
  Xc <- s * Xc
  
  # Add group column and column names
  group <- ifelse(x == 1, "A", "B")
  colnames(Xc) <- paste0("X", 1:ncol(Xc))
  Xc <- data.frame(Group = group, Xc)
  
  return(Xc)
}


# Replace 0 for transformation that can not handle 0
replace_zeros <- function(df, constant = 0.5) {
  df[df == 0] <- constant
  return(df)
}
# Additive log Ratio
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
# Additive Logit Ratio
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

# Additive Arcsine Ratio
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

# Centered Log Ratio
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

# Centered Logit Ratio
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

# Centered Arcsine Ratio
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

# Isometric Log Ratio
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

# Additive Dual Group Boxcox Ratio
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

# Centered Dual Group Boxcox Ratio
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
  data[, -1] <- data[, -1] / rowSums(data[, -1])
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

test_transformations <- function(simulated_data, transformation_functions) {
  results_list <- list()
  
  for (beta_s in names(simulated_data)) {
    
    data <- simulated_data[[beta_s]]
    results_list[[beta_s]] <- list()
    
    for (trans_name in names(transformation_functions)) {
      
      trans_func <- transformation_functions[[trans_name]]
      transformed_data <- trans_func(data) 
      group1_data <- transformed_data[transformed_data$Group == "A", -1]
      group2_data <- transformed_data[transformed_data$Group == "B", -1]
      
      p_values <- numeric(ncol(group1_data))
      
      for (i in 1:ncol(group1_data)) {
        compare_data <- data.frame(
          group = rep(c("A", "B"), each = nrow(group1_data)),
          datas = c(group1_data[, i], group2_data[, i])
        )
        test_result <- tryCatch({
          t.test(datas ~ group, data = compare_data)$p.value
        }, error = function(e) {
          # Return NA if there is an error
          return(NA)
        })
        if (is.na(test_result)) {
          p_values[i] <- NA
        } else {
          p_values[i] <- test_result
        }
      }
      
      # Adjust the p-values using the Bonferroni method
      adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
      
      results_list[[beta_s]][[trans_name]] <- adjusted_p_values
    }
  }
  return(results_list)
}

calculate_power_fdr <- function(results, ground_truth_results, significance_level = 0.05) {
  metrics_results <- list()
  
  for (beta_s in names(results)) {
    metrics_results[[beta_s]] <- list()
    for (trans_name in names(results[[beta_s]])) {
      p_values <- results[[beta_s]][[trans_name]]
      true_effects <- ground_truth_results[[beta_s]] 
      significant_count = sum(true_effects)
      TP <- sum(p_values < significance_level & true_effects)
      FP <- sum(p_values < significance_level & !true_effects)
      
      FN <- sum(p_values >= significance_level & true_effects)
      
      power <- TP / significant_count
      if (is.na(TP) || is.na(FP) || (TP + FP == 0)) {
        FDR <- NA  
      } else {
        FDR <- FP / (TP + FP)
      }
      
      metrics_results[[beta_s]][[trans_name]] <- list(Power = power, FDR = FDR)
    }
  }
  return(metrics_results)
}
generate_results_list <- function(num_datasets = 100, num_elements = 50) {
  results_list <- list()
  
  for (i in 1:num_datasets) {
    data_name <- paste0("data_", i)
    results_list[[data_name]] <- c(rep(TRUE, num_elements / 2), rep(FALSE, num_elements / 2))
  }
  
  return(results_list)
}
results_list <- generate_results_list()




options(warn=-1)
simulate_multiple_datasets <- function(num_datasets = 100, n = 100, p = 50, beta0 = 5, beta = 8, alpha = 3, q = 0.3) {
  datasets <- list()
  
  for (i in 1:num_datasets) {
    data_name <- paste0("data_", i)
    datasets[[data_name]] <- generate_data(n, p, beta0, beta, alpha, q)
  }
  
  return(datasets)
}

main_simulation <- function(num_datasets, n, p, beta0, beta, alpha, q) {
  simulated_data <- simulate_multiple_datasets(num_datasets, n, p, beta0, beta, alpha, q)
  transformation_functions <- list(
    "ALR" = alr_transformation,
    "CLR" = clr_transformation,
    "Additive Logit Ratio" = logit_alr_transformation,
    "Additive Arcsine Ratio" = arcsine_alr_transformation,
    "Centered Logit Ratio" = logit_clr_transformation,
    "Centered Arcsine Ratio" = arcsine_clr_transformation,
    "Additive (New)Power Ratio" = dual_group_boxcox_alr_transformation,
    "Centered (New)Power Ratio" = dual_group_boxcox_clr_transformation,
    "Boxcox in Ratio - William" = boxcox_alr_transformation_william,
    "New Boxcox in Ratio" = boxcox_alr_transformation_dual_group,
    "Additive New Logit Ratio" = new_logit_alr_transformation,
    "Centered New Logit Ratio" = new_logit_clr_transformation,
    "ILR" = ilr_transformation,
    "Additive Power Ratio" = power_alr_transformation,
    "Centered Power Ratio" = power_clr_transformation
  )
  
  results <- test_transformations(simulated_data, transformation_functions)
  ground_truth_results <- results_list # calculate_ground_truth(simulated_data)
  metrics_results <- calculate_power_fdr(results, ground_truth_results)
  summary_table <- tibble(beta_s = character(), transformation = character(), power = double(), fdr = double())
  
  for (beta_s in names(metrics_results)) {
    for (trans_name in names(metrics_results[[beta_s]])) {
      power <- metrics_results[[beta_s]][[trans_name]]$Power
      fdr <- metrics_results[[beta_s]][[trans_name]]$FDR
      
      summary_table <- rbind(summary_table, tibble(beta_s = beta_s, transformation = trans_name, power = power, fdr = fdr))
    }
  }
  
  return(summary_table)
}

run_multiple_simulations <- function(num_datasets, n, p) {
  all_results <- tibble()
  total_iterations <- 5 * 5 * 5 * length(seq(0, 0.9, by = 0.2))
  current_iteration <- 0
  for (alpha in seq(1, 10, by = 2)) {
    for (beta0 in seq(1, 10, by = 2)) {
      for (beta in seq(1, 10, by = 2)) {
        for (q in seq(0, 0.9, by = 0.2)) {
          result <- main_simulation(num_datasets, n, p, beta0, beta, alpha, q)
          result <- result %>%
            mutate(alpha = alpha, beta0 = beta0, beta = beta, q = q)
          all_results <- bind_rows(all_results, result)
          current_iteration <- current_iteration + 1
          cat("Progress:", current_iteration, "/", total_iterations, "\n")
        }
      }
    }
  }
  
  return(all_results)
}

# Parameters
num_datasets <- 100
n <- 100
p <- 50

# Run simulations
all_simulation_results <- run_multiple_simulations(num_datasets, n, p)

# Summarize results
summary_results <- all_simulation_results %>%
  group_by(alpha, beta0, beta, q, transformation) %>%
  summarize(
    mean_power = mean(power, na.rm = TRUE),
    sd_power = sd(power, na.rm = TRUE),
    mean_fdr = mean(fdr, na.rm = TRUE),
    sd_fdr = sd(fdr, na.rm = TRUE)
  )

write.csv(summary_results, "/home/yxz3116/summary_non_significant_as_reference.csv", row.names = FALSE)

