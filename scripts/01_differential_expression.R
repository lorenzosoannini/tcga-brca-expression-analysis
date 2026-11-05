# Historical exam workflow, with portable paths and explicit input checks.
# Run from the repository root. See docs/METHODOLOGY.md before interpreting FDR.
suppressPackageStartupMessages({library(stringr); library(pheatmap)})
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[1] else "data/matrix_RNAseq_BRCA.txt"
out <- if (length(args) >= 2) args[2] else "runs/BRCA"
if (!file.exists(input)) stop("Input missing. See data/README.md: ", input)
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory must be empty: ", out)
}
x <- as.matrix(read.delim(input, row.names = 1, check.names = FALSE, quote = ""))
if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0)) {
  stop("Expected finite, non-negative numeric expression values.")
}
if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) stop("Duplicate identifiers")
# Retain the first sample per patient and condition, as in the exam script.
select_samples <- function(pattern) {
  ids <- grep(pattern, colnames(x), value = TRUE)
  ids[!duplicated(str_extract(ids, "TCGA-\\w+-\\w+"))]
}
normal <- select_samples("TCGA-\\w+-\\w+-1\\d")
cancer <- select_samples("TCGA-\\w+-\\w+-0\\d")
patient <- function(z) str_extract(z, "TCGA-\\w+-\\w+")
common <- intersect(patient(normal), patient(cancer))
if (length(common) < 3) stop("At least three matched patient pairs are required.")
normal <- normal[match(common, patient(normal))]
cancer <- cancer[match(common, patient(cancer))]
stopifnot(identical(patient(normal), patient(cancer)))
x <- x[, c(normal, cancer), drop = FALSE]
x <- x[rowMeans(x) != 0, , drop = FALSE]
x <- log2(x + 1)
iqr <- apply(x, 1, IQR)
x <- x[iqr > quantile(iqr, 0.4), , drop = FALSE]
logFC <- rowMeans(x[, cancer, drop = FALSE]) - rowMeans(x[, normal, drop = FALSE])
# Deliberately preserve the historical fold-change prefilter and BH test family.
keep <- abs(logFC) >= 1
x <- x[keep, , drop = FALSE]; logFC <- logFC[keep]
if (!nrow(x)) stop("No genes passed the fold-change filter.")
pvalue <- vapply(seq_len(nrow(x)), function(i) {
  t.test(x[i, normal], x[i, cancer], paired = TRUE)$p.value
}, numeric(1))
padj <- p.adjust(pvalue, method = "BH")
keep <- is.finite(padj) & padj <= 0.01
x <- x[keep, , drop = FALSE]; logFC <- logFC[keep]
pvalue <- pvalue[keep]; padj <- padj[keep]
if (nrow(x) < 2) stop("Fewer than two genes passed the adjusted p-value filter.")
ids <- str_split_fixed(rownames(x), "\\|", 2)
deg <- data.frame(genes = ids[, 1], ensemblID = ids[, 2], pvalue = pvalue,
                  pval_adj = padj, logFC = logFC,
                  direction = ifelse(logFC < 0, "DOWN", "UP"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.table(deg, file.path(out, "DEG.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(x, file.path(out, "matrix_DEG.txt"), sep = "\t", col.names = NA, quote = FALSE)
for (kind in c("normal", "cancer")) {
  writeLines(get(kind), file.path(out, paste0("list_", kind, ".txt")))
}
pdf(file.path(out, "volcano_plot_pval.pdf"))
plot(logFC, -log10(padj), pch = 20, col = ifelse(logFC < 0, "#138A72", "#7C3EAD"),
     xlab = "Mean log2(expression + 1) difference", ylab = "-log10(BH-adjusted p-value)",
     main = "Selected genes: historical exam workflow")
abline(v = c(-1, 1), h = 2, lty = 2, col = "grey")
dev.off()
annotation <- data.frame(condition = rep(c("Normal", "Tumor"), each = length(common)),
                         row.names = colnames(x))
pheatmap(x, scale = "row", border_color = NA, show_rownames = FALSE, show_colnames = FALSE,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         clustering_method = "average", annotation_col = annotation,
         filename = file.path(out, "heatmap.pdf"), width = 10, height = 10)
writeLines(c(paste("Input:", normalizePath(input)), paste("Matched pairs:", length(common)),
             paste("Selected genes:", nrow(deg)), capture.output(sessionInfo())),
           file.path(out, "session-info.txt"))
message("Selected ", nrow(deg), " genes from ", length(common), " matched pairs. Output: ", out)
