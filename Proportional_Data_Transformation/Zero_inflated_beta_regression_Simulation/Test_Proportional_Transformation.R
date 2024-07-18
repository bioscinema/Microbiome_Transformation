# Define Transformation
library(MASS)
library(reticulate)
library(betareg)
library(moments)
library(dplyr)
source_python("Scource_Code.py")

log_transformation = function(x){
  return(log(x))
}
logit_transformation = function(x){
  return(log(x/(1-x)))
}
tangent_transformation = function(x){
  return(tan(pi*(x-0.5)))
}
arcsine_transformation = function(x){
  return(asin(sqrt(x)))
}

replace_zeros <- function(values, replacement_value = 1e-10) {
  values[values == 0] <- replacement_value
  return(values)
}
boxcox_transformation <- function(dataframe) {
  x = replace_zeros(dataframe[["value"]])
  boxcox_result <- boxcox(lm(x ~ 1, y = TRUE), plotit = FALSE)
  lambda_optimal <- boxcox_result$x[which.max(boxcox_result$y)]
  if (lambda_optimal == 0) {
    transformed_data <- log(x)
  } else {
    transformed_data <- (x^lambda_optimal - 1) / lambda_optimal
  }
  return(data.frame(group = dataframe$group, value = transformed_data))
}



Dual_Group_boxcox_transformation <- function(data, group_factor = 1) {
  factors = as.vector(unique(data[[group_factor]]))
  n = nrow(data)
  group1 = replace_zeros(data[data[group_factor] == factors[1], -group_factor, drop = FALSE])
  group2 = replace_zeros(data[data[group_factor] == factors[2], -group_factor, drop = FALSE])
  transformed = matrix(nrow = nrow(data), ncol = ncol(data) - 1)
  
  for(i in 1:(ncol(data) - 1)) {
    
    py$p1 = as.numeric(group1[, i])
    py$p2 = as.numeric(group2[, i])
    
    transformed_results <- py$dual_group_boxcox_transformation2(py$p1, py$p2)
    transformed[1:nrow(group1), i] <- transformed_results[[1]]
    transformed[(nrow(group1) + 1):n, i] <- transformed_results[[2]]
  }
  
  group_column = c(rep(factors[1], nrow(group1)), rep(factors[2], nrow(group2)))
  transformed_df <- data.frame(Group = group_column, transformed)
  colnames(transformed_df) = c("group", "value")
  return(transformed_df)
}

Dual_Group_logit_transformation <- function(data, group_factor = 1) {
  factors = as.vector(unique(data[[group_factor]]))
  n = nrow(data)
  group1 = data[data[group_factor] == factors[1], -group_factor, drop = FALSE]
  group2 = data[data[group_factor] == factors[2], -group_factor, drop = FALSE]
  transformed = matrix(nrow = nrow(data), ncol = ncol(data) - 1)
  
  for(i in 1:(ncol(data) - 1)) {
    
      py$p1 = as.numeric(group1[, i])
      py$p2 = as.numeric(group2[, i])
    
    transformed_results <- py$new_logit_transformation_s(py$p1, py$p2)
    transformed[1:nrow(group1), i] <- transformed_results[[1]]
    transformed[(nrow(group1) + 1):n, i] <- transformed_results[[2]]
  }
  
  group_column = c(rep(factors[1], nrow(group1)), rep(factors[2], nrow(group2)))
  transformed_df <- data.frame(Group = group_column, transformed)
  colnames(transformed_df) = c("group", "value")
  return(transformed_df)
}


generate_zero_inflated_beta_data <- function(n, p, eta_significant, eta_not_significant, nu_a, nu_b, x_a, x_b, beta, q) {
  # Calculate the number of observations per group
  n1 <- n / 2
  n2 <- n / 2
  
  # Create a grouping variable
  group <- factor(rep(c('A', 'B'), each = n1))
  
  # Initialize an empty data frame to store the results
  data <- data.frame(group = group)
  
  # Generate values for significant variables
  for (j in 1:(p/2)) {
    # Calculate mu for group A and B
    mu_a <- exp(eta_significant * x_a + beta) / (1 + exp(eta_significant * x_a + beta))
    mu_b <- exp(eta_significant * x_b + beta) / (1 + exp(eta_significant * x_b + beta))
    
    # Simulate data for group A
    group_a <- rbeta(n1, mu_a * nu_a, (1 - mu_a) * nu_a)
    
    # Simulate data for group B
    group_b <- rbeta(n2, mu_b * nu_b, (1 - mu_b) * nu_b)
    
    # Combine and apply zero-inflation
    combined <- c(group_a, group_b)
    combined <- combined * (rbinom(n, 1, 1 - q) == 1)
    
    # Bind the values to the main data frame
    data <- cbind(data, setNames(data.frame(combined), paste0('X', j)))
  }
  
  # Generate values for non-significant variables
  for (j in (p/2+1):p) {
    # Calculate mu for group A and B
    mu_a <- exp(eta_not_significant * x_a + beta) / (1 + exp(eta_not_significant * x_a + beta))
    mu_b <- exp(eta_not_significant * x_b + beta) / (1 + exp(eta_not_significant * x_b + beta))
    
    # Simulate data for group A
    group_a <- rbeta(n1, mu_a * nu_a, (1 - mu_a) * nu_a)
    
    # Simulate data for group B
    group_b <- rbeta(n2, mu_b * nu_b, (1 - mu_b) * nu_b)
    
    # Combine and apply zero-inflation
    combined <- c(group_a, group_b)
    combined <- combined * (rbinom(n, 1, 1 - q) == 1)
    
    # Bind the values to the main data frame
    data <- cbind(data, setNames(data.frame(combined), paste0('X', j)))
  }
  
  return(data)
}


simulate_data <- function(n, p, eta_significant, eta_not_significant, nu_a, nu_b, x_a, x_b, beta, q) {
  data = generate_zero_inflated_beta_data(n, p, eta_significant, eta_not_significant, nu_a, nu_b, x_a, x_b, beta, q)
  
  return(data)
}

perform_tests <- function(data, significance_level = 0.05, adjustment_method = "BH") {
  results <- list()
  p_values <- list()
  
  # T-test
  t_p_value <- t.test(value ~ group, data = data)$p.value
  p_values$t <- t_p_value
  
  # Logit-transformed T-test
  logit_data <- transform(data, value = logit_transformation(replace_zeros(value)))
  logit_t_p_value <- t.test(value ~ group, data = logit_data, na.action = na.exclude)$p.value
  p_values$logit_t <- logit_t_p_value
  
  # Log-transformed T-test
  log_data <- transform(data, value = log_transformation(replace_zeros(value)))
  log_t_p_value <- t.test(value ~ group, data = log_data, na.action = na.exclude)$p.value
  p_values$log_t <- log_t_p_value
  
  # Tangent-transformed T-test
  tangent_data <- transform(data, value = tangent_transformation(value))
  tangent_t_p_value <- t.test(value ~ group, data = tangent_data, na.action = na.exclude)$p.value
  p_values$tangent_t <- tangent_t_p_value
  
  # Arcsine-transformed T-test
  arcsine_data <- transform(data, value = arcsine_transformation(value))
  arcsine_t_p_value <- t.test(value ~ group, data = arcsine_data, na.action = na.exclude)$p.value
  p_values$arcsine_t <- arcsine_t_p_value
  
  # Wilcoxon Test
  wilcoxon_p_value <- wilcox.test(value ~ group, data = data)$p.value
  p_values$wilcoxon <- wilcoxon_p_value
  
  
  # Dual_group_logit-transformed T-test
  dual_group_logit_data <- Dual_Group_logit_transformation(data)
  dual_group_logit_t_p_value <- t.test(value ~ group, data = dual_group_logit_data, na.action = na.exclude)$p.value
  p_values$dual_group_logit_t <- dual_group_logit_t_p_value
  
  # Dual_group_boxcox-transformed T-test
  dual_group_boxcox_data <- Dual_Group_boxcox_transformation(data)
  dual_group_boxcox_t_p_value <- t.test(value ~ group, data = dual_group_boxcox_data, na.action = na.exclude)$p.value
  p_values$dual_group_boxcox_t <- dual_group_boxcox_t_p_value
  
  # Boxcox_transformation T-test
  boxcox_data <- boxcox_transformation(data)
  boxcox_t_p_value <- t.test(value ~ group, data = boxcox_data, na.action = na.exclude)$p.value
  p_values$boxcox_t <- boxcox_t_p_value
  
  
  # Adjust p-values
  adjusted_p_values <- p.adjust(unlist(p_values), method = adjustment_method)
  adjusted_p_values <- as.list(adjusted_p_values)
  
  # Determine significance based on adjusted p-values
  results$t <- adjusted_p_values$t < significance_level
  results$logit_t <- adjusted_p_values$logit_t < significance_level
  results$log_t <- adjusted_p_values$log_t < significance_level
  results$tangent_t <- adjusted_p_values$tangent_t < significance_level
  results$arcsine_t <- adjusted_p_values$arcsine_t < significance_level
  results$wilcoxon <- adjusted_p_values$wilcoxon < significance_level
  results$dual_group_logit_t <- adjusted_p_values$dual_group_logit_t < significance_level
  results$dual_group_boxcox_t <- adjusted_p_values$dual_group_boxcox_t < significance_level
  results$boxcox_t <- adjusted_p_values$boxcox_t < significance_level
  
  
  return(results)
}

# Function to simulate data and compute power and FDR for multiple tests

simulate_and_calculate_power_fdr <- function(n, p, eta_significant, eta_not_significant, nu_a, nu_b, x_a, x_b, beta, q,  num_simulations) {
  significance_level <- 0.05
  
  # Initialize lists to collect results for each test type
  results_list <- list()
  test_types <- c("t", "logit_t", "log_t", "tangent_t", "arcsine_t", "wilcoxon", "dual_group_logit_t","dual_group_boxcox_t",  "boxcox_t")
  
  # Initialize counts for true positives, false positives
  for (test_type in test_types) {
    results_list[[test_type]] <- list(
      true_positives = numeric(num_simulations),
      false_positives = numeric(num_simulations)
    )
  }
  
  # Simulation loop
  for (i in 1:num_simulations) {
    data <- simulate_data(n, p, eta_significant, eta_not_significant, nu_a, nu_b, x_a, x_b, beta, q)
    
    for (j in 1:p) {
      column_data <- data.frame(group = data$group, value = data[[paste0("X", j)]])
      colnames(column_data) <- c("group", "value")
      results <- perform_tests(column_data, significance_level)
      
      # Determine if the column should be significant (first 25 are significant)
      ground_truth_significant <- (j <= 25)
      
      # Update results for each test type
      for (test_type in test_types) {
        test_result <- results[[test_type]]
        
        if (isTRUE(test_result)) {
          if (ground_truth_significant) {
            results_list[[test_type]]$true_positives[i] <- results_list[[test_type]]$true_positives[i] + 1
          } else {
            results_list[[test_type]]$false_positives[i] <- results_list[[test_type]]$false_positives[i] + 1
          }
        }
      }
    }
  }
  
  # Calculate power, FDR, and their SDs
  power_results <- list()
  fdr_results <- list()
  power_sd <- list()
  fdr_sd <- list()
  
  for (test_type in test_types) {
    total_true_tests <- 25 * num_simulations
    total_false_tests <- 25 * num_simulations
    
    # Power is the ratio of true positives out of the total true tests
    true_positive_rates <- results_list[[test_type]]$true_positives / 25
    power_results[[test_type]] <- mean(true_positive_rates)
    power_sd[[test_type]] <- sd(true_positive_rates)
    
    # FDR is the ratio of false positives out of the total false tests
    false_positive_rates <- results_list[[test_type]]$false_positives / 25
    fdr_results[[test_type]] <- mean(false_positive_rates)
    fdr_sd[[test_type]] <- sd(false_positive_rates)
  }
  
  # Return the results
  return(list(
    power = power_results,
    fdr = fdr_results,
    power_sd = power_sd,
    fdr_sd = fdr_sd
  ))
}

simulate_and_combine_results <- function(n, p, eta_values, q_values, eta_not_significant, nu_a, nu_b, x_a, x_b, beta, num_simulations) {
  combined_results <- data.frame()
  
  for (eta in eta_values) {
    for (q in q_values) {
      result <- simulate_and_calculate_power_fdr(n, p, eta, eta_not_significant, nu_a, nu_b, x_a, x_b, beta, q, num_simulations)
      df <- data.frame(
        Method = names(result$power),
        Power = unlist(result$power),
        FDR = unlist(result$fdr),
        Power_sd = unlist(result$power_sd),
        FDR_sd = unlist(result$fdr_sd),
        eta = eta,
        q = q
      )
      combined_results <- rbind(combined_results, df)
    }
  }
  
  return(combined_results)
}

n <- 100
p <- 50
eta_values <- c(-0.7, -0.5)
q_values <- seq(0, 0.7, by = 0.1)
eta_not_significant <- 0
nu_a <- 5
nu_b <- 5
x_a <- 1
x_b <- 2
beta <- -2
num_simulations <- 100

# Generate combined results
combined_results <- simulate_and_combine_results(n, p, eta_values, q_values, eta_not_significant, nu_a, nu_b, x_a, x_b, beta, num_simulations)

# Save to CSV
write.csv(combined_results, "Power_summary.csv", row.names = FALSE)




adjust_values <- function(values) {
  epsilon <- 1e-5
  pmax(epsilon, pmin(values, 1 - epsilon))
}
calculate_skewness_kurtosis <- function(data) {
  results <- list()
  
  original_stats <- sapply(data[-1], function(col) {
    list(
      skewness = list(
        A = skewness(col[data$group == "A"], na.rm = TRUE),
        B = skewness(col[data$group == "B"], na.rm = TRUE)
      ),
      kurtosis = list(
        A = kurtosis(col[data$group == "A"], na.rm = TRUE),
        B = kurtosis(col[data$group == "B"], na.rm = TRUE)
      )
    )
  })
  
  results$original <- original_stats
  
  transformations_no_group <- list(
    logit = logit_transformation,
    log = log_transformation,
    tangent = tangent_transformation,
    arcsine = arcsine_transformation
  )
  
  for (transformation_name in names(transformations_no_group)) {
    transformed_data <- data.frame(group = data$group)
    
    for (colname in colnames(data)[-1]) {
      transformed_data[[colname]] <- transformations_no_group[[transformation_name]](adjust_values(data[[colname]]))
    }
    
    transformed_stats <- sapply(transformed_data[-1], function(col) {
      list(
        skewness = list(
          A = skewness(col[transformed_data$group == "A"], na.rm = TRUE),
          B = skewness(col[transformed_data$group == "B"], na.rm = TRUE)
        ),
        kurtosis = list(
          A = kurtosis(col[transformed_data$group == "A"], na.rm = TRUE),
          B = kurtosis(col[transformed_data$group == "B"], na.rm = TRUE)
        )
      )
    })
    
    results[[transformation_name]] <- transformed_stats
  }
  
  transformations_with_group <- list(
    dual_group_logit = Dual_Group_logit_transformation,
    dual_group_boxcox = Dual_Group_boxcox_transformation,
    boxcox = boxcox_transformation
  )
  
  for (transformation_name in names(transformations_with_group)) {
    transformed_data <- data.frame(group = data$group)
    
    for (colname in colnames(data)[-1]) {
      transformed_column_data <- data.frame(group = data$group, value = adjust_values(data[[colname]]))
      transformed_column_data <- transformations_with_group[[transformation_name]](transformed_column_data)
      transformed_data[[colname]] <- transformed_column_data$value
    }
    
    transformed_stats <- sapply(transformed_data[-1], function(col) {
      list(
        skewness = list(
          A = skewness(col[transformed_data$group == "A"], na.rm = TRUE),
          B = skewness(col[transformed_data$group == "B"], na.rm = TRUE)
        ),
        kurtosis = list(
          A = kurtosis(col[transformed_data$group == "A"], na.rm = TRUE),
          B = kurtosis(col[transformed_data$group == "B"], na.rm = TRUE)
        )
      )
    })
    
    results[[transformation_name]] <- transformed_stats
  }
  
  return(results)
}

# Function to summarize statistics
summarize_statistics <- function(values) {
  c(
    mean = mean(values, na.rm = TRUE),
    sd = sd(values, na.rm = TRUE)
  )
}

# Main simulation function
perform_simulations <- function(eta, q, n_simulations) {
  all_stats <- vector("list", n_simulations)
  
  for (i in 1:n_simulations) {
    data <- generate_zero_inflated_beta_data(n, p, eta, 0, nu_a, nu_b, x_a, x_b, beta, q)
    all_stats[[i]] <- calculate_skewness_kurtosis(data)
  }
  
  transformation_names <- names(all_stats[[1]])
  data_names <- c("skewness", "kurtosis")
  group_names <- c("A", "B")
  values <- list()
  
  for (i in seq_len(n_simulations)) {
    for (transformation_name in transformation_names) {
      for (data_name in data_names) {
        for (group_name in group_names) {
          temp <- c()
          for (j in seq_len(ncol(all_stats[[i]][[transformation_name]]))) {
            subdata <- all_stats[[i]][[transformation_name]][, j]
            temp <- c(temp, subdata[[data_name]][[group_name]])
          }
          if (is.null(values[[transformation_name]])) values[[transformation_name]] <- list()
          if (is.null(values[[transformation_name]][[data_name]])) values[[transformation_name]][[data_name]] <- list()
          if (is.null(values[[transformation_name]][[data_name]][[group_name]])) values[[transformation_name]][[data_name]][[group_name]] <- c()
          values[[transformation_name]][[data_name]][[group_name]] <- c(values[[transformation_name]][[data_name]][[group_name]], mean(temp))
        }
      }
    }
  }
  
  summarize_table <- list()
  for (transformation_name in transformation_names) {
    summarize_table[[transformation_name]] <- list()
    for (data_name in data_names) {
      summarize_table[[transformation_name]][[data_name]] <- list()
      for (group_name in group_names) {
        summarize_table[[transformation_name]][[data_name]][[group_name]] <- summarize_statistics(values[[transformation_name]][[data_name]][[group_name]])
      }
    }
  }
  
  return(summarize_table)
}

# Perform simulations for each eta and q combination
results <- list()

eta_values <- c(-0.5, -0.7)
q_values <- c(0, 0.3, 0.5, 0.7)
n_simulations <- 100

for (eta in eta_values) {
  for (q in q_values) {
    simulation_result <- perform_simulations(eta, q, n_simulations)
    results[[paste0("eta_", eta, "_q_", q)]] <- simulation_result
  }
}

# Convert results to data frame for CSV output
df_list <- list()

for (result_name in names(results)) {
  result <- results[[result_name]]
  
  # Extract eta and q values correctly
  eta_value <- as.numeric(gsub("eta_([^_]+)_q_.*", "\\1", result_name))
  q_value <- as.numeric(gsub(".*_q_([^_]+)", "\\1", result_name))
  
  for (transformation_name in names(result)) {
    for (data_name in names(result[[transformation_name]])) {
      for (group_name in names(result[[transformation_name]][[data_name]])) {
        temp <- result[[transformation_name]][[data_name]][[group_name]]
        df_list[[length(df_list) + 1]] <- data.frame(
          eta = eta_value,
          q = q_value,
          transformation = transformation_name,
          metric = data_name,
          group = group_name,
          mean = temp["mean"],
          sd = temp["sd"]
        )
      }
    }
  }
}
final_df <- bind_rows(df_list)
write.csv(final_df, "Skewness_kurtosis.csv", row.names = FALSE)
