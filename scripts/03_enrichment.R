# Optional live Enrichr queries. Database revisions may differ from the exam exports.
suppressPackageStartupMessages(library(enrichR))
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[1] else "results/exam/BRCA/DEG.txt"
out <- if (length(args) >= 2) args[2] else "runs/enrichment"
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory must be empty: ", out)
}
deg <- read.delim(input, check.names = FALSE)
stopifnot(all(c("genes", "direction") %in% names(deg)))
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021", "KEGG_2021_Human", "DisGeNET")
available <- listEnrichrDbs()
if (!all(dbs %in% available$libraryName)) stop("An exam database is unavailable; use archived exports.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (direction in c("UP", "DOWN")) {
  genes <- unique(deg$genes[deg$direction == direction])
  if (!length(genes)) next
  tables <- enrichr(genes, dbs)
  for (db in dbs) {
    tab <- tables[[db]]
    if (is.null(tab) || !"Adjusted.P.value" %in% names(tab)) stop("Missing enrichment response: ", db)
    write.table(tab, file.path(out, paste0(direction, "_", db, "_all.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    tab <- tab[is.finite(tab$Adjusted.P.value) & tab$Adjusted.P.value < 0.05, , drop = FALSE]
    tab <- tab[order(tab$Adjusted.P.value), , drop = FALSE]
    write.table(tab, file.path(out, paste0(direction, "_", db, "_significant.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
writeLines(c(paste("Query time (UTC):", format(Sys.time(), tz = "UTC")), capture.output(sessionInfo())),
           file.path(out, "session-info.txt"))
