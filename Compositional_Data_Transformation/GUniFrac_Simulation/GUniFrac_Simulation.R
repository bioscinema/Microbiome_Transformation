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
packages <- c("compositions", "reticulate", "MASS", "dplyr", "tibble", "GUniFrac")

# Install and load each package
for (pkg in packages) {
  install_and_load(pkg, user_lib)
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


rdirichlet.m <- function (alpha) {
  Gam <- matrix(rgamma(length(alpha), shape = alpha), nrow(alpha), ncol(alpha))
  t(t(Gam) / colSums(Gam))
}
EstPara <- function (ref.otu.tab) {
  
  if (is.null(rownames(ref.otu.tab))) {
    rownames(ref.otu.tab) <- paste0('OTU', 1 : nrow(ref.otu.tab))
  } # otu * sample
  samplenames = colnames(ref.otu.tab)
  taxnames = rownames(ref.otu.tab)
  
  dirmult.paras <- dirmult::dirmult(t(ref.otu.tab))
  
  gamma = dirmult.paras$gamma
  names(gamma) = names(dirmult.paras$pi)
  
  # Add pseduo count(each OTU add gamma estimated from dirmult)
  ref.otu.tab = sapply(1:ncol(ref.otu.tab), function (i) gamma + ref.otu.tab[,i]) # C_ij otu * sample
  
  # back to Dirichlet, calculate the true proportion
  ref.otu.tab.p <- rdirichlet.m(ref.otu.tab) # P_ij nOTU*nSam
  colnames(ref.otu.tab.p) = samplenames
  rownames(ref.otu.tab.p) = taxnames
  
  # order OTUs by mean OTU proportion, for later selection
  ord = order(rowMeans(ref.otu.tab.p), decreasing = TRUE)
  ref.otu.tab.p =  ref.otu.tab.p[ord,]
  
  # apply size factor
  Si = exp(rnorm(ncol(ref.otu.tab.p)))
  ref.otu.tab0 = t(t(ref.otu.tab.p)*Si)
  colnames(ref.otu.tab0) = colnames(ref.otu.tab.p)
  return(list(mu = ref.otu.tab.p, ref.otu.tab = ref.otu.tab0))
}
SimulateMSeqU<-function (para = EstPara, nSam = 100, nOTU = 500, diff.otu.pct = 0.1, 
                         diff.otu.direct = c("balanced", "unbalanced"), 
                         diff.otu.mode = c("abundant",  "rare", "mix","user_specified"), 
                         user_specified_otu = NULL,
                         covariate.type = c("binary", "continuous"), 
                         grp.ratio = 1, covariate.eff.mean = 1, covariate.eff.sd = 0, 
                         confounder.type = c("none", "binary", "continuous", "both"), 
                         conf.cov.cor = 0.6, conf.diff.otu.pct = 0, conf.nondiff.otu.pct = 0.1, 
                         confounder.eff.mean = 0, confounder.eff.sd = 0, error.sd = 0, 
                         depth.mu = 10000, depth.theta = 5, depth.conf.factor = 0, cont.conf, epsilon)
{
  ## Estimated Dirichlet Parameters
  model.paras <- para
  
  ## Select number of OTUs
  ref.otu.tab <- model.paras$ref.otu.tab[(1:(nOTU)), ]
  
  ## Select number of Samples
  idx.otu <- rownames(ref.otu.tab)
  idx.sample <- colnames(model.paras$ref.otu.tab)[1:nSam]
  ref.otu.tab = ref.otu.tab[, idx.sample]
  
  ## Confounding variable
  if (confounder.type == "none") {
    confounder.type <- "continuous"
    confounder.eff.mean <- 0
    confounder.eff.sd <- 0
    Z <- NULL
  }
  if (confounder.type == "continuous") {
    Z <- cbind(cont.conf)
  }
  
  if (confounder.type == "binary") {
    Z <- cbind(c(rep(0, nSam%/%2), rep(1, nSam - nSam%/%2)))
  }
  
  if (confounder.type == "both") {
    Z <- cbind(rnorm(nSam), c(rep(0, nSam%/%2), rep(1, nSam - nSam%/%2)))
  } 
  
  ## Covariate of Interest
  rho <- sqrt(conf.cov.cor^2/(1 - conf.cov.cor^2))
  
  if (covariate.type == "continuous") {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + epsilon#rnorm(nSam)
  }
  
  if (covariate.type == "binary") {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + epsilon#rnorm(nSam)
    X <- cbind(ifelse(X <= quantile(X, grp.ratio/(1 + grp.ratio)), 0, 1))
  }
  
  rownames(X) <- colnames(ref.otu.tab)
  covariate.eff.mean1 = covariate.eff.mean
  covariate.eff.mean2 = covariate.eff.mean
  
  ## Simulate OTU Absolute Abundance
  ## Balanced OTU Counts
  if (diff.otu.direct == "balanced") {
    if (diff.otu.mode == "user_specified") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                           rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*% 
        t(scale(X))
    }
    if (diff.otu.mode == "abundant") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                           rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*% 
        t(scale(X))
    }
    else if (diff.otu.mode == "rare") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                           rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*% 
        t(scale(X))
    }
    else {
      eta.diff <- c(sample(c(rnorm(floor(nOTU/4), mean = -covariate.eff.mean1, sd = covariate.eff.sd), 
                             rnorm(floor(nOTU/2) - floor(nOTU/4), mean = covariate.eff.mean1, sd = covariate.eff.sd))), 
                    sample(c(rnorm(floor((nOTU - floor(nOTU/2))/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                             rnorm(nOTU - floor(nOTU/2) - floor((nOTU - floor(nOTU/2))/2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*% 
        t(scale(X))
    }
  }
  
  ## Unbalanced OTU Counts
  if (diff.otu.direct == "unbalanced") {
    if (diff.otu.mode == "user_specified") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*% 
        t(scale(X))
    }
    if (diff.otu.mode == "abundant") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*% 
        t(scale(X))
    }
    else if (diff.otu.mode == "rare") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*% 
        t(scale(X))
    }
    else {
      eta.diff <- c(sample(c(rnorm(floor(nOTU/2), mean = covariate.eff.mean1, sd = covariate.eff.sd))), 
                    sample(c(rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*%
        t(scale(X))
    }
  }
  eta.conf <- sample(c(rnorm(floor(nOTU/2), mean = -confounder.eff.mean, sd = confounder.eff.sd), 
                       rnorm(nOTU - floor(nOTU/2), mean = confounder.eff.mean, sd = confounder.eff.sd))) %*% 
    t(scale(scale(Z) %*% rep(1, ncol(Z))))
  
  otu.ord <- 1:(nOTU)
  diff.otu.ind <- NULL
  diff.otu.num <- round(diff.otu.pct * nOTU)
  
  if (diff.otu.mode == "user_specified") 
    diff.otu.ind <- c(diff.otu.ind, which(idx.otu %in% user_specified_otu))
  if (diff.otu.mode == "mix") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord, diff.otu.num))
  if (diff.otu.mode == "abundant") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[1:round(length(otu.ord)/4)], diff.otu.num))
  if (diff.otu.mode == "rare") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[round(3 * length(otu.ord)/4):length(otu.ord)], diff.otu.num))
  
  if (length(diff.otu.ind) >= round(nOTU * conf.diff.otu.pct)) {
    conf.otu.ind1 <- sample(diff.otu.ind, round(nOTU * conf.diff.otu.pct))
  } else {
    conf.otu.ind1 <- diff.otu.ind
  }
  conf.otu.ind <- c(conf.otu.ind1, 
                    sample(setdiff(1:(nOTU), diff.otu.ind), round(conf.nondiff.otu.pct * nOTU)))
  
  ## Calculate the new composition
  eta.diff[setdiff(1:(nOTU), diff.otu.ind), ] <- 0
  eta.conf[setdiff(1:(nOTU), conf.otu.ind), ] <- 0
  eta.error <- matrix(rnorm(nOTU * nSam, 0, error.sd), nOTU, nSam)
  eta.exp <- exp(t(eta.diff + eta.conf + eta.error))
  eta.exp <- eta.exp * t(ref.otu.tab)
  ref.otu.tab.prop <- eta.exp/rowSums(eta.exp)
  ref.otu.tab.prop <- t(ref.otu.tab.prop)
  
  ## Simulate Sequencing Depth using Negative Binomial Distribution
  nSeq <- rnegbin(nSam, 
                  mu = depth.mu * exp(scale(X) * depth.conf.factor), 
                  theta = depth.theta)
  
  ## Simulate Absolute Abundance of OTU
  otu.tab.sim <- sapply(1:ncol(ref.otu.tab.prop), 
                        function(i) rmultinom(1, nSeq[i], ref.otu.tab.prop[, i]))
  colnames(otu.tab.sim) <- rownames(eta.exp)
  rownames(otu.tab.sim) <- rownames(ref.otu.tab)
  diff.otu.ind = (1:nOTU) %in% diff.otu.ind
  conf.otu.ind = (1:nOTU) %in% conf.otu.ind
  
  ## Output
  return(list(otu.tab.sim = otu.tab.sim, covariate = X, confounder = Z, 
              diff.otu.ind = diff.otu.ind, otu.names = idx.otu, conf.otu.ind = conf.otu.ind))
}

load("para1.RData")


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
      if (is.na(TP) || is.na(FP) || (TP + FP) == 0) {
        FDR <- NA
      } else {
        FDR <- FP / (TP + FP)
      }
      
      metrics_results[[beta_s]][[trans_name]] <- list(Power = power, FDR = FDR)
    }
  }
  
  return(metrics_results)
}


main <- function() {
  simulated_data <- gunifrac_list
  transformation_functions <- list(
    "ALR" = alr_transformation,
    "CLR" = clr_transformation,
    "Additive Logit Ratio" =logit_alr_transformation,
    "Additive Arcsine Ratio" = arcsine_alr_transformation,
    "Centered Logit Ratio" = logit_clr_transformation,
    "Centered Arcsine Ratio" = arcsine_clr_transformation,
    "Additive (New)Power Ratio" =dual_group_boxcox_alr_transformation,
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
  ground_truth_results <- results_list #calculate_ground_truth(simulated_data)  
  metrics_results <- calculate_power_fdr(results, ground_truth_results)  
  summary_table <- tibble(beta_s = character(), transformation = character(), power_fdr = character())
  
  for (beta_s in names(metrics_results)) {
    for (trans_name in names(metrics_results[[beta_s]])) {
      power <- metrics_results[[beta_s]][[trans_name]]$Power
      fdr <- metrics_results[[beta_s]][[trans_name]]$FDR
      
      summary_table <- rbind(summary_table, tibble(beta_s = beta_s, transformation = trans_name, power = power, fdr = fdr))
    }
  }
  
  return(summary_table)
}


library(dplyr)
options(warn=-1)
# Define the parameters to vary in the simulations
n_simulation = 100
param_combinations <- expand.grid(
  diff_otu_direct = c("unbalanced", "balanced"),
  diff_otu_mode = c("rare", "mix", "abundant"),
  depth_mu = c(10000, 1000, 100, 10),
  depth_theta = c(5, 10, 15),
  covariate_eff_sd = c(0, 0.5),
  confounder_eff_sd = c(0, 0.5),
  depth_conf_factor = c(0, 0.5)
)

# Initialize list to store all summarized results
all_summarized_results_list = list()
total_param_sets <- nrow(param_combinations)
# Run simulations with different parameters
for (param_set in seq(nrow(param_combinations))) {
  diff_otu_direct <- param_combinations$diff_otu_direct[param_set]
  diff_otu_mode <- param_combinations$diff_otu_mode[param_set]
  depth_mu <- param_combinations$depth_mu[param_set]
  depth_theta <- param_combinations$depth_theta[param_set]
  covariate_eff_sd <- param_combinations$covariate_eff_sd[param_set]
  confounder_eff_sd <- param_combinations$confounder_eff_sd[param_set]
  depth_conf_factor <- param_combinations$depth_conf_factor[param_set]
  
  gunifrac_list = list()
  results_list = list()
  sim_info_list = list()
  
  for (i in 1:n_simulation) {
    Simulated_data <- SimulateMSeqU(
      para = para1, nSam = 100, nOTU = 100, diff.otu.pct = 0.1,
      diff.otu.direct = diff_otu_direct,
      diff.otu.mode = diff_otu_mode,
      user_specified_otu = NULL,
      covariate.type = "binary",
      grp.ratio = 1, covariate.eff.mean = 1, covariate.eff.sd = covariate_eff_sd,
      confounder.type = "both",
      conf.cov.cor = 0.6, conf.diff.otu.pct = 0, conf.nondiff.otu.pct = 0.1,
      confounder.eff.mean = 0, confounder.eff.sd = confounder_eff_sd, error.sd = 0,
      depth.mu = depth_mu, depth.theta = depth_theta, depth.conf.factor = depth_conf_factor, cont.conf = 0, epsilon = 0
    )
    
    # Extract the OTU table and metadata
    otu_df <- as.data.frame(Simulated_data$otu.tab.sim)
    diff_otu_ids <- Simulated_data$otu.names[Simulated_data$diff.otu.ind]
    meta.dat <- data.frame(
      X = Simulated_data$covariate, 
      Z1 = Simulated_data$confounder[, 1],
      Z2 = Simulated_data$confounder[, 2]
    )
    meta.dat$sample <- rownames(meta.dat)
    
    # Transpose the OTU table
    t_otu_df <- as.data.frame(t(otu_df))
    
    # Keep the count data and merge with metadata
    t_otu_df$sample <- rownames(t_otu_df)
    merged_data_current <- merge(t_otu_df, meta.dat, by = "sample")
    merged_data_current$group <- ifelse(merged_data_current$X == 0, "A", "B")
    
    # Process the OTU table to include the group column and remove unnecessary columns
    otu_df_processed <- merged_data_current %>%
      select(group, everything()) %>%
      select(-c(X, Z1, Z2, sample))
    
    # Rename columns
    col_names <- colnames(otu_df_processed)[-1]
    indices <- match(diff_otu_ids, col_names)
    indications <- ifelse(1:100 %in% indices, TRUE, FALSE)
    results_list[[paste0(diff_otu_direct, "_", diff_otu_mode, "_data_", i)]] <- indications
    col_number <- as.integer(ncol(otu_df_processed))
    colnames(otu_df_processed) <- c("group", paste0("X", seq(1, col_number - 1, by = 1)))
    gunifrac_list[[paste0(diff_otu_direct, "_", diff_otu_mode, "_data_",i)]] = otu_df_processed
    
    # Add simulation info to the dataframe
    sim_info <- data.frame(
      simulation = paste0("sim_", i),
      diff_otu_direct = diff_otu_direct,
      diff_otu_mode = diff_otu_mode,
      depth_mu = depth_mu,
      depth_theta = depth_theta,
      covariate_eff_sd = covariate_eff_sd,
      confounder_eff_sd = confounder_eff_sd,
      depth_conf_factor = depth_conf_factor
    )
    sim_info_list[[paste0("info_", diff_otu_direct, "_", diff_otu_mode, "_", i)]] = sim_info
  }
  
  # Combine all simulation data into a single data frame for the current parameter set
  combined_data <- do.call(rbind, lapply(gunifrac_list, function(x) data.frame(x)))
  combined_info <- do.call(rbind, lapply(sim_info_list, function(x) data.frame(x)))
  final_data <- cbind(combined_info, combined_data)
  
  
  
  # Run the main() function to get the result for the current parameter set
  simulated_result <- main()
  
  # Summarize the result
  summarized_result <- simulated_result %>%
    group_by(transformation) %>%
    summarise(
      mean_power = mean(power, na.rm = TRUE),
      sd_power = sd(power, na.rm = TRUE),
      mean_fdr = mean(fdr, na.rm = TRUE),
      sd_fdr = sd(fdr, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Add the parameters as the first columns of the summarized result
  summarized_result_with_params <- cbind(
    data.frame(
      diff_otu_direct = diff_otu_direct,
      diff_otu_mode = diff_otu_mode,
      depth_mu = depth_mu,
      depth_theta = depth_theta,
      covariate_eff_sd = covariate_eff_sd,
      confounder_eff_sd = confounder_eff_sd,
      depth_conf_factor = depth_conf_factor,
      stringsAsFactors = FALSE
    ),
    summarized_result
  )
  
  # Store the summarized result in the results list
  all_summarized_results_list[[paste0(diff_otu_direct, "_", diff_otu_mode, "_", depth_mu, "_", depth_theta, "_", covariate_eff_sd, "_", confounder_eff_sd, "_", depth_conf_factor)]] <- summarized_result_with_params
  
  print(paste("Completed parameter set", param_set, "of", total_param_sets))
}

# Combine all summarized results into a single data frame
all_summarized_results_df <- do.call(rbind, all_summarized_results_list)

# Save all summarized results to a CSV file, excluding the beta_s column if present
final_results <- all_summarized_results_df %>%
  select(-contains("beta_s"))

write.csv(final_results, "/home/yxz3116/all_summarized_results_GUniFrac.csv", row.names = FALSE)
