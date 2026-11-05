rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("~/Desktop/Bioinf")
dirRes <- "results/"
if (!file.exists(dirRes)) {
  dir.create(dirRes)
}
#############################################
library(stringr)
library(pheatmap) 
#############################################
tumor <- "BRCA"
dirTumor <- paste0(dirRes, tumor, "/")
if (!file.exists(dirTumor)) {
  dir.create(dirTumor)
} 
#############################################
# input data
filename_in <- paste0("matrix_RNAseq_", tumor, ".txt")

# output data
filename_list_normal <- paste0(dirTumor, "list_normal.txt")
filename_list_cancer <- paste0(dirTumor, "list_cancer.txt") 
filename_DEG <- paste0(dirTumor, "DEG.txt")
filename_matrix_DEG <- paste0(dirTumor, "matrix_DEG.txt")
filename_heatmap <- paste0(dirTumor, "heatmap.pdf")
#############################################
# parameters
prc_IQR <- 0.4
thr_fc <- 2
thr_pval <- 0.01
paired <- TRUE
#############################################
# import and read data
tmp <- read.table(filename_in, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, nrows = 5)
classes <- sapply(tmp, class)
tmp <- read.table(filename_in, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "",colClasses = classes)

genes <- rownames(tmp)
samples <- colnames(tmp)

# filter matched patients in case and control groups
pzN <- grep('TCGA-\\w+-\\w+-1\\d', samples, value = TRUE) # normal samples
pzC <- grep('TCGA-\\w+-\\w+-0\\d', samples, value = TRUE) # cancer samples

pzN <- pzN[!duplicated(str_extract(pzN, 'TCGA-\\w+-\\w+'))]
pzC <- pzC[!duplicated(str_extract(pzC, 'TCGA-\\w+-\\w+'))]

nameN <- str_extract(pzN, 'TCGA-\\w+-\\w+')
nameC <- str_extract(pzC, 'TCGA-\\w+-\\w+')

common_pz <- intersect(nameN, nameC) # common patients

pzN_com <- unlist(lapply(common_pz, function(x){grep(x, pzN, value = TRUE)})) # common normal samples barcodes
pzC_com <- unlist(lapply(common_pz, function(x){grep(x, pzC, value = TRUE)})) # common cancer samples barcodes

dataN <- tmp[, pzN_com] # common normal samples data
dataC <- tmp[, pzC_com] # common cancer samples data

data <- cbind(dataN, dataC)
pz_com <- colnames(data)

# cleanup
rm(tmp, samples, pzN, pzC, nameN, nameC, common_pz, classes)
#############################################
# PRE_PROCESSING
#############################################
# remove genes not expressed both in normal and cancer samples
overall_mean <- apply(data, 1, mean)
ind <- which(overall_mean == 0)
dataN <- dataN[-ind, ]
dataC <- dataC[-ind, ]
data <- data[-ind, ]
genes <- genes[-ind]

# logarithmic transformation
dataN <- log2(dataN + 1)
dataC <- log2(dataC + 1) 
data <- log2(data + 1)

# compute IQR
variation <- apply(data, 1, IQR)

# IQR filtering
thr_prc <- quantile(variation, prc_IQR)
ind <- which(variation <= thr_prc)
dataN <- dataN[-ind, ]
dataC <- dataC[-ind, ]
data <- data[-ind, ]
genes <- genes[-ind]
# cleanup
rm(ind)

# plot IQR distribution
pdf(paste0(dirTumor, "IQR_distribution.pdf"))
hist(variation, breaks = 100, col = "lightblue", main = "IQR frequency distribution", xlab = "IQR value", ylab = "Frequency")
abline(v = thr_prc, col = "red", lwd = 2, lty = 2)
dev.off()

# compute and plot IQR after filtering
variation <- apply(data, 1, IQR)
pdf(paste0(dirTumor, "IQR_distribution_filtered.pdf"))
hist(variation, breaks = 100, col = "lightblue", main = "IQR frequency distribution", xlab = "IQR value", ylab = "Frequency")
dev.off()
#############################################
# FILTERING
#############################################
# compute fold change
pdf(paste0(dirTumor, "logFC_distribution.pdf"))
logFC <- rowMeans(dataC) - rowMeans(dataN)
hist(logFC, breaks = 100,col = "orange", main = "Fold change frequency distribution", xlab = "log FC", ylab = "Frequency")
abline(v = c(-log2(thr_fc), log2(thr_fc)), col = "grey", lwd = 3, lty = 2)
dev.off()

# FC filtering
ind <- which(abs(logFC) < log2(thr_fc)) 
if (length(ind) > 0) {
  dataN <- dataN[-ind, ]
  dataC <- dataC[-ind, ]
  data <- data[-ind, ]
  genes <- genes[-ind]
  logFC <- logFC[-ind]
}
# cleanup
rm(ind)

# compute p-value with paired t-test (#samples > 30 so we can assume normal distribution)
N <- ncol(dataN)
M <- ncol(dataC)
pval <- apply(data, 1, function(x){t.test(x[1:N], x[(N+1):(N+M)], paired = paired)$p.value})
# compute adjusted p-value
pval_adj <- p.adjust(pval, method = "fdr")

# p-value filtering
ind <- which(pval_adj > thr_pval)

if (length(ind) > 0) {
  dataN <- dataN[-ind, ]
  dataC <- dataC[-ind, ]
  data <- data[-ind, ]
  genes <- genes[-ind]
  logFC <- logFC[-ind]
  pval <- pval[-ind]
  pval_adj <- pval_adj[-ind]
}
# cleanup
rm(ind)

# volcano plot p-value distribution, green points on the left are down-regulated genes, purple points on the right are up-regulated genes
pdf(paste0(dirTumor, "volcano_plot_pval.pdf"))
plot(logFC, -log10(pval_adj), pch = 20, col = ifelse(logFC < 0, "green", "purple"), main = "Volcano plot", xlab = "log FC", ylab = "-log10 p-value")
abline(v = c(-log2(thr_fc), log2(thr_fc)), col = "grey", lwd = 3, lty = 2)
abline(h = -log10(thr_pval), col = "red", lwd = 3, lty = 2)
dev.off()

# plot the up/down-regulated genes frequency with a pie chart
pdf(paste0(dirTumor, "pie_chart.pdf"))
pie(table(ifelse(logFC < 0, "DOWN", "UP")), col = c("green", "purple"), main = "Up/down-regulated genes")
# add percentages
pct <- round(table(ifelse(logFC < 0, "DOWN", "UP"))/length(logFC) * 100)
lbls <- paste(names(pct), pct)
lbls <- paste(lbls, "%", sep = "")
legend("topright", legend = lbls, cex = 0.8, fill = c("green", "purple"))
dev.off()

#############################################
# EXPORT RESULTS
#############################################
direction <- ifelse(logFC < 0, "DOWN", "UP")

DEG <- data.frame(str_split_fixed(genes, "\\|", 2))
colnames(DEG) <- c("geneSymbol", "ensemblID")

results <- data.frame(genes = DEG$geneSymbol, ensemblID = DEG$ensemblID, pvalue = pval, pval_adj = pval_adj, logFC = logFC, direction = direction)

write.table(results, file = filename_DEG, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data, file = filename_matrix_DEG, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(pzN_com, file = filename_list_normal, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(pzC_com, file = filename_list_cancer, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#############################################
# OTHER PLOTS
#############################################



#############################################
# ADD BOXPLOT OF THE MOST INTERESTING GENES!!
#############################################



# boxplot the most up and down-regulated genes
# (notch represents the 95% confidence interval around the median)
pdf(paste0(dirTumor, "boxplot_max.pdf"))
ind <- which.max(results$logFC)
gene_id <- results$genes[ind]

df <- data.frame(normal = t(dataN[ind, ]), cancer = t(dataC[ind, ]), row.names = NULL)
colnames(df) <- c(paste(gene_id, "normal"), paste(gene_id, "cancer"))
boxplot(df, main = paste0(gene_id, ", ", "adjusted p-value = ", format(pval_adj[ind], digits = 2)), notch = TRUE,
        ylab = "Gene expression value", xlab = "Condition", col = c("green", "red"), pars = list(boxwex = 0.3, staplewex = 0.6))

dev.off() 

pdf(paste0(dirTumor, "boxplot_min.pdf"))
ind <- which.min(results$logFC)
gene_id <- results$genes[ind]

df <- data.frame(normal = t(dataN[ind, ]), cancer = t(dataC[ind, ]), row.names = NULL)
colnames(df) <- c(paste(gene_id, "normal"), paste(gene_id, "cancer"))
boxplot(df, main = paste0(gene_id, ", ", "adjusted p-value = ", format(pval_adj[ind], digits = 2)), notch = TRUE,
        ylab = "Gene expression value", xlab = "Condition", col = c("green", "red"), pars = list(boxwex = 0.3, staplewex = 0.6))

dev.off()
#cleanup
rm(ind, gene_id, df)
#############################################
# HEATMAP
#############################################
test <- grepl('TCGA-\\w+-\\w+-1\\d', colnames(data))
condition <- ifelse(test, "normal", "cancer")

annotation <- data.frame(condition = condition, row.names = colnames(data))
vect_color <- c("green", "red")
names(vect_color) <- unique(condition)
annotation_colors <- list(condition = vect_color)

pheatmap(data, scale = "row", border_color = NA, cluster_cols = T, cluster_rows = T,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", clustering_method = "average",
         annotation_col = annotation, annotation_colors = annotation_colors, color = colorRampPalette(colors = c("blue", "blue3", "black", "yellow", "yellow3"))(100),
         show_rownames = F, show_colnames = F, cutree_cols = 2, cutree_rows = 2, width = 10, height = 10, filename = filename_heatmap)

