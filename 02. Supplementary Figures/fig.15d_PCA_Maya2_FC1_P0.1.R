library(tidyverse)
library(ggbiplot)
library(ggplot2)

# Load the data
Total_Maya2_root <- read.csv('Total_Maya2_root.csv')

# Select relevant columns
PCA_table <- Total_Maya2_root %>%
  select(2, 5:20)

# List up each gene
WT_list <- PCA_table %>%
  filter(abs(Log2FC_WT) > 1 & Padj_WT < 0.1) %>%
  pull(Gene)

fls2_list <- PCA_table %>%
  filter(abs(Log2FC_fls2) > 1 & Padj_fls2 < 0.1) %>%
  pull(Gene)

# Make unique gene list
uni_list <- union(WT_list, fls2_list)

# Generate PCA table
PCA_table_f <- PCA_table %>%
  filter(Gene %in% uni_list)

# Select relevant columns for PCA
PCA_table_f_s <- PCA_table_f %>%
  select(1, 6:17)

# Add unique sample identifiers to avoid duplicate row names
PCA_table_f_s$SampleID <- make.unique(as.character(PCA_table_f_s$Gene))

# Set SampleID as rownames and remove the Gene column
PCA_table_f_s_unique <- PCA_table_f_s %>%
  column_to_rownames(var = "SampleID") %>%
  select(-Gene)

# Transpose the data frame
PCA_table_f_s_t <- PCA_table_f_s_unique %>%
  t() %>%
  as.data.frame()

# Add group information to rownames
group_info <- c(rep("DMSO_WT", 3), rep("Maya1_WT", 3), rep("DMSO_fls2", 3), rep("Maya1_fls2", 3))

# Perform PCA
pca_result <- prcomp(PCA_table_f_s_t, center = TRUE, scale. = TRUE)

# Check PCA result
summary(pca_result)

# Define groups as factor
groups <- factor(group_info)

# Plot PCA
ggbiplot(pca_result, groups = groups, var.axes = FALSE, ellipse = TRUE) +
  geom_point(aes(color = groups, shape = groups), size = 3) +
  scale_color_manual(values = c('DMSO_WT'='#92c5de', 'Maya1_WT'='#0571b0', 'DMSO_fls2'='#f4a582', 'Maya1_fls2'='#ca0020')) +
  scale_shape_manual(values = c('DMSO_WT' = 15, 'Maya1_WT' = 16, 'DMSO_fls2' = 17, 'Maya1_fls2' = 18)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.ticks = element_line(colour = "black"),
        plot.margin = unit(c(1, 1, 1, 1), "line"),
        axis.title = element_text(size = 0),
        legend.position = "none",
        legend.title = element_text(size = 0),
        legend.spacing.x = unit(0.5, "cm"),
        legend.key.size = unit(0.9, "cm"),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 0))


