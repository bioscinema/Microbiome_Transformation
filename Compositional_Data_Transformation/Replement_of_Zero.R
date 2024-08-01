library(dplyr)

generate_data <- function(n = 100, p = 500, beta0 = 3, beta = 3, alpha = 2, q = 0.3) {
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
  
  pi <- replicate(p, rbinom(n, 1, 1 - q))
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

replace_zeros <- function(df, constant = 0.5) {
  df[df == 0] <- constant
  return(df)
}

alr_transformation <- function(data, group_factor = 1, ref_component = ncol(data)-1, constant = 0.5) {
  data <- replace_zeros(data, constant)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  factors <- as.vector(unique(data[[group_factor]]))
  group1 <- data[data[group_factor] == factors[1], -group_factor]
  group2 <- data[data[group_factor] == factors[2], -group_factor]
  group1_alr <- log(group1[, -ref_component] / group1[, ref_component])
  group2_alr <- log(group2[, -ref_component] / group2[, ref_component])
  colnames(group1_alr) <- colnames(data[, -c(group_factor, (ref_component+1))])
  colnames(group2_alr) <- colnames(data[, -c(group_factor, (ref_component+1))])
  group1_alr <- cbind(Group = factors[1], group1_alr)
  group2_alr <- cbind(Group = factors[2], group2_alr)
  combined <- rbind(group1_alr, group2_alr)
  return(combined)
}

clr_transformation <- function(data, group_factor = 1, constant = 0.5) {
  data <- replace_zeros(data, constant)
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  factors <- as.vector(unique(data[[group_factor]]))
  group1 <- data[data[group_factor] == factors[1], -group_factor]
  group2 <- data[data[group_factor] == factors[2], -group_factor]
  gm_group1 <- exp(rowMeans(log(group1)))
  gm_group2 <- exp(rowMeans(log(group2)))
  clr_group1 <- log(group1 / gm_group1)
  clr_group2 <- log(group2 / gm_group2)
  clr_group1 <- cbind(Group = factors[1], clr_group1)
  clr_group2 <- cbind(Group = factors[2], clr_group2)
  combined <- rbind(clr_group1, clr_group2)
  return(combined)
}

# Function to perform t-tests and calculate power and FDR
analyze_transformed_data <- function(transformed_data, n_significant = 250) {
  p_values <- sapply(transformed_data[, -1], function(x) {
    t.test(x ~ transformed_data$Group)$p.value
  })
  
  # Calculate power and FDR
  significant <- p_values[1:n_significant] < 0.05
  non_significant <- p_values[(n_significant+1):length(p_values)] < 0.05
  power <- sum(significant) / n_significant
  fdr <- sum(non_significant) / (sum(non_significant) + sum(significant))
  
  return(data.frame(Power = power, FDR = fdr))
}

# Run simulations
set.seed(123)
n_simulations <- 300
q_values <- c(0, 0.3, 0.5, 0.7)
constant_values <- c(0.1, 0.5, 1, 2,5)
results_list <- list()

for (q in q_values) {
  for (constant in constant_values) {
    alr_results_list <- list()
    clr_results_list <- list()
    
    for (i in 1:n_simulations) {
      data <- generate_data(q = q)
      
      alr_data <- alr_transformation(data, constant = constant)
      clr_data <- clr_transformation(data, constant = constant)
      
      alr_results <- analyze_transformed_data(alr_data)
      clr_results <- analyze_transformed_data(clr_data)
      
      alr_results_list[[i]] <- alr_results
      clr_results_list[[i]] <- clr_results
    }
    
    # Calculate mean power and FDR
    alr_results_mean <- do.call(rbind, alr_results_list) %>%
      summarise(across(everything(), mean))
    clr_results_mean <- do.call(rbind, clr_results_list) %>%
      summarise(across(everything(), mean))
    
    # Combine results into a table
    results_table <- rbind(data.frame(Method = "ALR", Q = q, Constant = constant, alr_results_mean),
                           data.frame(Method = "CLR", Q = q, Constant = constant, clr_results_mean))
    
    results_list[[paste0("Q_", q, "_Constant_", constant)]] <- results_table
  }
}

final_results <- do.call(rbind, results_list)
# ALR
alr_data <- final_results %>% filter(Method == "ALR")
alr_power_anova <- aov(Power ~ Q + Constant, data = alr_data)
alr_fdr_anova <- aov(FDR ~ Q + Constant, data = alr_data)

# CLR
clr_data <- final_results %>% filter(Method == "CLR")
clr_power_anova <- aov(Power ~ Q + Constant, data = clr_data)
clr_fdr_anova <- aov(FDR ~ Q + Constant, data = clr_data)


summary(alr_power_anova)
summary(alr_fdr_anova)
summary(clr_power_anova)
summary(clr_fdr_anova)