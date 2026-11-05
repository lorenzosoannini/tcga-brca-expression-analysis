# PCA on already-selected genes: descriptive, not independent validation.
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[1] else "results/exam/BRCA"
out <- if (length(args) >= 2) args[2] else "runs/PCA"
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory must be empty: ", out)
}
x <- as.matrix(read.delim(file.path(input, "matrix_DEG.txt"), row.names = 1,
                          check.names = FALSE, quote = ""))
normal <- readLines(file.path(input, "list_normal.txt"))
cancer <- readLines(file.path(input, "list_cancer.txt"))
samples <- c(cancer, normal)
stopifnot(!anyDuplicated(samples), all(samples %in% colnames(x)), all(is.finite(x)))
pca <- prcomp(t(x[, samples, drop = FALSE]), center = TRUE, scale. = TRUE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
scores <- data.frame(sample = samples, PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     condition = c(rep("Tumor", length(cancer)), rep("Normal", length(normal))))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
g <- ggplot(scores, aes(PC1, PC2, color = condition)) + geom_point(alpha = 0.7) +
  stat_ellipse(level = 0.68, type = "norm") + theme_bw() +
  labs(x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2]),
       color = "Condition", title = "PCA on genes selected in the exam analysis")
ggsave(file.path(out, "score_plot.pdf"), g, width = 6, height = 5)
write.table(scores, file.path(out, "scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(component = seq_along(variance), explained_percent = variance),
            file.path(out, "explained_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pca$rotation, file.path(out, "loadings.tsv"), sep = "\t", quote = FALSE, col.names = NA)
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
message(sprintf("PC1 = %.2f%%; PC2 = %.2f%%. Output: %s", variance[1], variance[2], out))
