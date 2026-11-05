# Install packages for differential expression and PCA; add optional enrichment.
packages <- c("stringr", "pheatmap", "ggplot2")
if ("--enrichment" %in% commandArgs(trailingOnly = TRUE)) packages <- c(packages, "enrichR")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
