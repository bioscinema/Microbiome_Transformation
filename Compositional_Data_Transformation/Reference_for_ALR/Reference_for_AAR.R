library(phyloseq) 
library(DESeq2)
library(dplyr)
library(tibble) 
library(tidyverse)
library(tidyverse)
library(ggplot2)



load("physeq.RData")


filter_rows <- function(data) {
  filtered_data <- data[apply(data, 1, function(row) {
    sum(row != 0) >= 5
  }), ]
  return(filtered_data)
}





# Get significant by using deseq2
otu_data <- as.data.frame(otu_table(physeq.tree))
filtered_otu_data <- otu_data[, grepl("\\.MD", colnames(otu_data))]
filtered_otu_data <- filter_rows(filtered_otu_data)
sample_data <- as.data.frame(sample_data(physeq.tree))
filtered_sample_data <- sample_data[grepl("\\.MD", rownames(sample_data)), ]
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(filtered_otu_data),
  colData = filtered_sample_data,
  design = ~ group
)
dds <- DESeq(dds)
res <- results(dds)
significant_results <- res %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  filter(padj < 0.05)
deseq2_significant <- as.vector(significant_results["OTU"])$OTU






arcsine <- function(x) {
  return(asin(sqrt(x)))
}
arcsine_alr_transformation <- function(data, group_factor = 1, ref_component) {
  data[, -1] <- data[, -1] / rowSums(data[, -1])
  data_without_group = as.matrix(data[,-group_factor])
  p = dim(data_without_group)[2]
  T = rbind(diag(p-1),rep(-1,p-1))
  temp <- T[ref_component, ]       # Store the first row in a temporary variable
  T[ref_component, ] <- T[p, ]  # Replace the first row with the last row
  T[p, ] <- temp
  t_data = arcsine(data_without_group) %*% T
  transformed_data = data.frame(data[,group_factor], t_data)
  colnames(transformed_data)[1] <- "Group"
  colnames(transformed_data)[-1] <- colnames(data_without_group[,-ncol(data_without_group)])
  return(transformed_data)
}

# Test significant for ALR Transformation
## Situation 1: replace 0 with 0.5, without remove outliers
otu_table_data <- otu_table(physeq.tree)
otu_table_transposed <- t(filter_rows(otu_table_data))
data = as.data.frame(otu_table_transposed)
filtered_df <- data %>%
  filter(grepl("\\.MD", rownames(data)))
df <- filtered_df %>%
  mutate(Group = case_when(
    grepl("^LTS", rownames(.)) ~ "A",
    grepl("^STS", rownames(.)) ~ "B"
  ))
df_0 <- df %>%
  select(Group, everything())
df_cleaned <- (df_0)


otu_table_compositional <- df_cleaned %>%
  select(-Group) %>%
  as.matrix()
otu_table_compositional <- (otu_table_compositional) / rowSums((otu_table_compositional),na.rm = TRUE)
otu_table_compositional <- as.data.frame(otu_table_compositional)
otu_table_compositional$Group <- df_cleaned$Group
otu_table_compositional <- otu_table_compositional %>%
  select(Group, everything())

identify_significant_columns <- function(data, group_factor = 1) {
  results <- data.frame(Reference = integer(), SignificantColumns = character(), stringsAsFactors = FALSE)
  num_cols <- ncol(data) - 1
  
  for (ref_component in 1:num_cols) {
    transformed_data <- arcsine_alr_transformation(data, group_factor, ref_component)
    
    group1_data <- transformed_data %>% filter(Group == unique(transformed_data$Group)[1])
    group2_data <- transformed_data %>% filter(Group == unique(transformed_data$Group)[2])
    
    p_values <- sapply(2:ncol(transformed_data), function(col_index) {
      if (col_index == ref_component + 1) {
        return(NA)
      } else {
        p_value <- tryCatch({
          t.test(group1_data[[col_index]], group2_data[[col_index]])$p.value
        }, error = function(e) {
          return(NA)
        })
      }
    })
    significant_cols <- which(p.adjust(p_values, method = "BH", n = length(p_values)) < 0.05)
    
    results <- rbind(results, data.frame(Reference = ref_component, SignificantColumns = paste(significant_cols, collapse = ", ")))
  }
  
  return(results)
}
df <- (otu_table_compositional)
significant_columns_results <- identify_significant_columns(df, group_factor = 1)


num_cols <- ncol(df) - 1
result_matrix <- matrix("white", nrow = num_cols, ncol = num_cols)
for (i in 1:nrow(significant_columns_results)) {
  ref <- significant_columns_results$Reference[i]
  sig_cols <- as.numeric(unlist(strsplit(significant_columns_results$SignificantColumns[i], ", ")))
  if (length(sig_cols) > 0 && !all(is.na(sig_cols))) {
    result_matrix[ref, sig_cols] <- "red"
  }
  result_matrix[ref, ref] <- "blue"
}
result_df <- as.data.frame(result_matrix)
result_df$Reference <- 1:num_cols
plot_df <- pivot_longer(result_df, cols = -Reference, names_to = "Column", values_to = "Color")
plot_df$Column <- as.numeric(gsub("V", "", plot_df$Column))
plot_df$Color <- factor(plot_df$Color, levels = c("white", "blue", "red"))
values_to_annotate <- which(colnames(df) %in% deseq2_significant) - 1
triangles_df_x <- data.frame(Column = values_to_annotate, Reference = rep(1, length(values_to_annotate)))
triangles_df_y <- data.frame(Column = rep(1, length(values_to_annotate)), Reference = values_to_annotate)
g <- ggplot(plot_df, aes(x = Column, y = Reference, fill = Color)) +
  geom_tile() +
  scale_fill_manual(values = c("white" = "white", "blue" = "blue", "red" = "red"),
                    labels = c("Non-Significant", "Reference", "Significant")) +
  geom_point(data = triangles_df_x, aes(x = Column, y = Reference), shape = 17, color = "green", size = 1, inherit.aes = FALSE) +
  geom_point(data = triangles_df_y, aes(x = Column, y = Reference), shape = 17, color = "green", size = 1, inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "Significant Columns by Reference Component for AAR Transformation",
       x = "Variables",
       y = "Variables",
       fill = "Legend") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), aspect.ratio = 1)
plot_width <- 10
plot_height <- 10
ggsave(g, file = "Arcsine_Different_Reference_for_outlier.eps", device = "eps", width = plot_width, height = plot_height)

