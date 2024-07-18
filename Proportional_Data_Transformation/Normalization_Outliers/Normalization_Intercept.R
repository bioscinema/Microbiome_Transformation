# Load necessary library
library(MASS)

# Define transformations as functions
logit <- function(x) {
  log(x / (1 - x))
}

arcsine_transform <- function(x) {
  (2 / pi) * asin(sqrt(x))
}

boxcox_transform <- function(x) {
  boxcox_result <- boxcox(x ~ 1, plot = FALSE)
  lambda <- boxcox_result$x[which.max(boxcox_result$y)]
  if (lambda == 0) {
    return(list(transformed = log(x), lambda = lambda))
  } else {
    return(list(transformed = (x^lambda - 1) / lambda, lambda = lambda))
  }
}

set.seed(123)

n <- 100
num_simulations <- 100000
beta0_values <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

# Initialize storage for results
results_no_outlier <- data.frame(
  beta0 = numeric(),
  estimated_intercept_orig = numeric(),
  se_orig = numeric(),
  estimated_intercept_log = numeric(),
  se_log = numeric(),
  estimated_intercept_logit = numeric(),
  se_logit = numeric(),
  estimated_intercept_arcsine = numeric(),
  se_arcsine = numeric(),
  estimated_intercept_boxcox = numeric(),
  se_boxcox = numeric()
)
lambda_values_no_outlier <- list()
# Perform simulations
for (beta0 in beta0_values) {
  intercept_orig <- numeric(num_simulations)
  intercept_log <- numeric(num_simulations)
  intercept_logit <- numeric(num_simulations)
  intercept_arcsine <- numeric(num_simulations)
  intercept_boxcox <- numeric(num_simulations)
  lambda_boxcox <- numeric(num_simulations)
  for (i in 1:num_simulations) {
    # Generate original data
    epsilon <- runif(n, -0.15, 0.15)
    y_orig <- beta0 + epsilon
    
    # Fit models
    fit_orig <- lm(y_orig ~ 1)
    fit_log <- lm(log(y_orig) ~ 1)
    fit_logit <- lm(logit(y_orig) ~ 1)
    fit_arcsine <- lm(arcsine_transform(y_orig) ~ 1)
    y_boxcox <- boxcox_transform(y_orig)
    fit_boxcox <- lm(y_boxcox$transformed ~ 1)
    
    # Store intercept estimates
    intercept_orig[i] <- coef(fit_orig)[1]
    intercept_log[i] <- coef(fit_log)[1]
    intercept_logit[i] <- coef(fit_logit)[1]
    intercept_arcsine[i] <- coef(fit_arcsine)[1]
    intercept_boxcox[i] <- coef(fit_boxcox)[1]
    lambda_boxcox[i] <- y_boxcox$lambda
  }
  
  # Compute mean and standard error of intercept estimates
  mean_orig <- mean(intercept_orig)
  se_orig <- sd(intercept_orig)
  mean_log <- mean(intercept_log)
  se_log <- sd(intercept_log)
  mean_logit <- mean(intercept_logit)
  se_logit <- sd(intercept_logit)
  mean_arcsine <- mean(intercept_arcsine)
  se_arcsine <- sd(intercept_arcsine)
  mean_boxcox <- mean(intercept_boxcox)
  se_boxcox <- sd(intercept_boxcox)
  
  # Append results to data frame
  results_no_outlier <- rbind(results_no_outlier, data.frame(
    beta0 = beta0,
    estimated_intercept_orig = mean_orig,
    se_orig = se_orig,
    estimated_intercept_log = mean_log,
    se_log = se_log,
    estimated_intercept_logit = mean_logit,
    se_logit = se_logit,
    estimated_intercept_arcsine = mean_arcsine,
    se_arcsine = se_arcsine,
    estimated_intercept_boxcox = mean_boxcox,
    se_boxcox = se_boxcox
  ))
  lambda_values_no_outlier[[as.character(beta0)]] <- lambda_boxcox
}

n <- 100 # Increased sample size due to outliers
num_simulations <- 100000
beta0_values <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
outliers <- c(0.9, 0.99, 0.999, 0.9999, 0.99999, 0.99999)

results_outlier <- data.frame( 
  beta0 = numeric(),
  estimated_intercept_orig = numeric(),
  se_orig = numeric(),
  estimated_intercept_log = numeric(),
  se_log = numeric(),
  estimated_intercept_logit = numeric(),
  se_logit = numeric(),
  estimated_intercept_arcsine = numeric(),
  se_arcsine = numeric(),
  estimated_intercept_boxcox = numeric(),
  se_boxcox = numeric()
)

# Perform simulations
for (beta0 in beta0_values) {
  intercept_orig <- numeric(num_simulations)
  intercept_log <- numeric(num_simulations)
  intercept_logit <- numeric(num_simulations)
  intercept_arcsine <- numeric(num_simulations)
  intercept_boxcox <- numeric(num_simulations)
  
  for (i in 1:num_simulations) {
    # Generate original data
    epsilon <- runif(n - length(outliers), -0.15, 0.15)
    y_orig <- c(beta0 + epsilon, outliers)
    
    # Fit models
    fit_orig <- lm(y_orig ~ 1)
    fit_log <- lm(log(y_orig) ~ 1)
    fit_logit <- lm(logit(y_orig) ~ 1)
    fit_arcsine <- lm(arcsine_transform(y_orig) ~ 1)
    y_boxcox <- boxcox_transform(y_orig)
    fit_boxcox <- lm(y_boxcox$transformed ~ 1)
    
    # Store intercept estimates
    intercept_orig[i] <- coef(fit_orig)[1]
    intercept_log[i] <- coef(fit_log)[1]
    intercept_logit[i] <- coef(fit_logit)[1]
    intercept_arcsine[i] <- coef(fit_arcsine)[1]
    intercept_boxcox[i] <- coef(fit_boxcox)[1]
  }
  
  # Compute mean and standard error of intercept estimates
  mean_orig <- mean(intercept_orig)
  se_orig <- sd(intercept_orig)
  mean_log <- mean(intercept_log)
  se_log <- sd(intercept_log)
  mean_logit <- mean(intercept_logit)
  se_logit <- sd(intercept_logit)
  mean_arcsine <- mean(intercept_arcsine)
  se_arcsine <- sd(intercept_arcsine)
  mean_boxcox <- mean(intercept_boxcox)
  se_boxcox <- sd(intercept_boxcox)
  
  # Append results to data frame
  results_outlier <- rbind(results_outlier, data.frame(
    beta0 = beta0,
    estimated_intercept_orig = mean_orig,
    se_orig = se_orig,
    estimated_intercept_log = mean_log,
    se_log = se_log,
    estimated_intercept_logit = mean_logit,
    se_logit = se_logit,
    estimated_intercept_arcsine = mean_arcsine,
    se_arcsine = se_arcsine,
    estimated_intercept_boxcox = mean_boxcox,
    se_boxcox = se_boxcox
  ))
}

library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

# Define a function to create the line plot with error bars and add horizontal lines
create_line_plot_with_errorbars <- function(df, y_var, y_label, se_var, title = NULL, show_legend = TRUE, y_limits = NULL, hlines = NULL, hlines_color = NULL) {
  df_melt <- melt(df, id.vars = "beta0", measure.vars = y_var)
  se_melt <- melt(df, id.vars = "beta0", measure.vars = se_var)
  df_melt$variable <- factor(df_melt$variable, levels = c("estimated_intercept_orig", "estimated_intercept_log", "estimated_intercept_logit", "estimated_intercept_arcsine", "estimated_intercept_boxcox"),
                             labels = c("Original", "Log", "Logit", "Arcsine", "Box-Cox"))
  se_melt$variable <- factor(se_melt$variable, levels = c("se_orig", "se_log", "se_logit", "se_arcsine", "se_boxcox"),
                             labels = c("Original", "Log", "Logit", "Arcsine", "Box-Cox"))
  df_combined <- merge(df_melt, se_melt, by = c("beta0", "variable"), suffixes = c("", "_se"))
  
  p <- ggplot(df_combined, aes(x = beta0, y = value, color = variable)) +
    geom_line() +
    geom_errorbar(aes(ymin = value - value_se, ymax = value + value_se), width = 0.02) +
    labs(x = expression(beta[0]), y = y_label) +
    theme_minimal() +
    ggtitle(title) +
    theme(legend.position = ifelse(show_legend, "right", "none"),
          plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(size = 8)) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5),
          axis.text = element_text(size = 9, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"))
  if (!is.null(y_limits)) {
    p <- p + ylim(y_limits)
  }
  if (!is.null(hlines)) {
    for (i in seq_along(hlines)) {
      p <- p + geom_hline(yintercept = hlines[[i]]$intercept, linetype = "dotted", color = hlines[[i]]$color)
    }
  }
  return(p)
}

# Determine y-axis limits for intercept
intercept_y_limits <- range(c(
  results_no_outlier$estimated_intercept_orig - results_no_outlier$se_orig, results_no_outlier$estimated_intercept_log - results_no_outlier$se_log,
  results_no_outlier$estimated_intercept_logit - results_no_outlier$se_logit, results_no_outlier$estimated_intercept_arcsine - results_no_outlier$se_arcsine, results_no_outlier$estimated_intercept_boxcox - results_no_outlier$se_boxcox,
  results_outlier$estimated_intercept_orig - results_outlier$se_orig, results_outlier$estimated_intercept_log - results_outlier$se_log,
  results_outlier$estimated_intercept_logit - results_outlier$se_logit, results_outlier$estimated_intercept_arcsine - results_outlier$se_arcsine, results_outlier$estimated_intercept_boxcox - results_outlier$se_boxcox,
  results_no_outlier$estimated_intercept_orig + results_no_outlier$se_orig, results_no_outlier$estimated_intercept_log + results_no_outlier$se_log,
  results_no_outlier$estimated_intercept_logit + results_no_outlier$se_logit, results_no_outlier$estimated_intercept_arcsine + results_no_outlier$se_arcsine, results_no_outlier$estimated_intercept_boxcox + results_no_outlier$se_boxcox,
  results_outlier$estimated_intercept_orig + results_outlier$se_orig, results_outlier$estimated_intercept_log + results_outlier$se_log,
  results_outlier$estimated_intercept_logit + results_outlier$se_logit, results_outlier$estimated_intercept_arcsine + results_outlier$se_arcsine, results_outlier$estimated_intercept_boxcox + results_outlier$se_boxcox
))

# Extract intercept values for beta0 = 0.2 and beta0 = 0.8 from no outlier data
hlines <- list(
  list(intercept = results_no_outlier$estimated_intercept_orig[results_no_outlier$beta0 == 0.2], color = "red"),
  list(intercept = results_no_outlier$estimated_intercept_orig[results_no_outlier$beta0 == 0.8], color = "red"),
  list(intercept = results_no_outlier$estimated_intercept_log[results_no_outlier$beta0 == 0.2], color = "darkgreen"),
  list(intercept = results_no_outlier$estimated_intercept_log[results_no_outlier$beta0 == 0.8], color = "darkgreen"),
  list(intercept = results_no_outlier$estimated_intercept_logit[results_no_outlier$beta0 == 0.2], color = "green"),
  list(intercept = results_no_outlier$estimated_intercept_logit[results_no_outlier$beta0 == 0.8], color = "green"),
  list(intercept = results_no_outlier$estimated_intercept_arcsine[results_no_outlier$beta0 == 0.2], color = "blue"),
  list(intercept = results_no_outlier$estimated_intercept_arcsine[results_no_outlier$beta0 == 0.8], color = "blue"),
  list(intercept = results_no_outlier$estimated_intercept_boxcox[results_no_outlier$beta0 == 0.2], color = "purple"),
  list(intercept = results_no_outlier$estimated_intercept_boxcox[results_no_outlier$beta0 == 0.8], color = "purple")
)

# Create plots for intercept with error bars and horizontal lines
plot_intercept_no_outlier <- create_line_plot_with_errorbars(results_no_outlier, c("estimated_intercept_orig", "estimated_intercept_log", "estimated_intercept_logit", "estimated_intercept_arcsine", "estimated_intercept_boxcox"), "Intercept", c("se_orig", "se_log", "se_logit", "se_arcsine", "se_boxcox"), "No Outliers", show_legend = TRUE, y_limits = intercept_y_limits, hlines = hlines)
plot_intercept_outlier <- create_line_plot_with_errorbars(results_outlier, c("estimated_intercept_orig", "estimated_intercept_log", "estimated_intercept_logit", "estimated_intercept_arcsine", "estimated_intercept_boxcox"), NULL, c("se_orig", "se_log", "se_logit", "se_arcsine", "se_boxcox"), "Outliers", show_legend = FALSE, y_limits = intercept_y_limits, hlines = hlines) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(5, 15, 5, 5))

# Extract legend from one of the plots
get_legend <- function(my_plot) {
  tmp <- ggplotGrob(my_plot + theme(legend.position = "right", legend.text = element_text(size = 8)))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(plot_intercept_no_outlier)

# Combine the plots and the legend in a grid layout
combined_plot <- arrangeGrob(
  arrangeGrob(plot_intercept_no_outlier + theme(legend.position = "none", plot.margin = margin(6, 0, 6, 6)), 
              plot_intercept_outlier + theme(legend.position = "none", plot.margin = margin(6, 6, 6, 0)), 
              nrow = 1, ncol = 2),
  legend,
  nrow = 1,
  top = textGrob("Intercept Estimates with Standard Error for Different Conditions", gp = gpar(fontsize = 11, fontface = "bold")),
  layout_matrix = rbind(c(1, 1, 2)),
  widths = c(1, 1, 0.3)
)

# Save the combined plot as EPS
ggsave("variance_and_outliers.eps", plot = combined_plot, device = "eps", width = 12, height = 5)
