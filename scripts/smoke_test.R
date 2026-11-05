# Exercise pairing, zero filtering, selection and output integrity on synthetic data.
set.seed(17)
workspace <- tempfile("brca-smoke-")
dir.create(workspace)
n <- 40; g <- 120
patients <- sprintf("TCGA-ZZ-%04d", seq_len(n))
normal <- paste0(patients, "-11A-01R-TEST-01")
cancer <- paste0(patients, "-01A-01R-TEST-01")
base <- matrix(runif(g * n, 4, 7), nrow = g)
shift <- rep(c(2, -2, 0), each = 40)
tumor <- base + shift + matrix(rnorm(g * n, sd = 0.2), nrow = g)
x <- cbind(2^base - 1, 2^tumor - 1)
rownames(x) <- paste0("GENE", seq_len(g), "|ENSG", seq_len(g))
colnames(x) <- c(normal, cancer)
x <- rbind(x, "ZERO|ENSGZERO" = rep(0, 2 * n))
# Deliberately shuffle columns so a successful run requires patient alignment.
x <- x[, sample(seq_len(ncol(x)))]
input <- file.path(workspace, "matrix.txt")
write.table(x, input, sep = "\t", quote = FALSE, col.names = NA)
output <- file.path(workspace, "expression")
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("scripts/01_differential_expression.R", shQuote(input), shQuote(output)))
stopifnot(status == 0)
deg <- read.delim(file.path(output, "DEG.txt"))
stopifnot(nrow(deg) > 20, all(deg$pval_adj <= 0.01), all(abs(deg$logFC) >= 1),
          all(c("UP", "DOWN") %in% deg$direction), !"ZERO" %in% deg$genes)
# The selected synthetic genes must follow the known direction of their shift.
idx <- as.integer(sub("GENE", "", deg$genes))
stopifnot(all(sign(deg$logFC) == sign(shift[idx])))
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("scripts/02_pca.R", shQuote(output), shQuote(file.path(workspace, "pca"))))
stopifnot(status == 0)
message("PASS: shuffled matched samples, synthetic effect directions, filtering, expression and PCA execution")
