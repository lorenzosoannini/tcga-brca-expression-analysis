# Optional: running the scripts

The exam results are already included in the repository. These commands are for
readers who want to explore the portable scripts; they are not required to review
the project and are not run by the repository's automated workflow.

Run from the repository root. The portable scripts were checked with R 4.4.3;
the historical package environment was not recorded. Dependencies are installed
from CRAN, rather than pinned to the original exam versions.

```bash
Rscript scripts/install_dependencies.R

# Explore PCA using the bundled selected-gene matrix.
Rscript scripts/02_pca.R

# Run differential expression only if you have the original input matrix.
Rscript scripts/01_differential_expression.R /path/to/matrix_RNAseq_BRCA.txt runs/BRCA
Rscript scripts/02_pca.R runs/BRCA runs/PCA-rerun
```

See [input requirements](../data/README.md) before supplying an expression matrix.
The original source-data download release and normalization unit were not fully
recorded; a current TCGA download should not be assumed to reproduce the exam.

Optional live enrichment:

```bash
Rscript scripts/install_dependencies.R --enrichment
Rscript scripts/03_enrichment.R results/exam/BRCA/DEG.txt runs/enrichment
```

This command needs internet access and submits gene symbols to Enrichr. Live
results may differ from archived exports because database contents can change.

Commands require an empty output directory, write under `runs/` by default, and
record their R session. They do not overwrite the archived exam outputs.
`smoke_test.R` is an optional execution check with synthetic data and known
expression shifts; it does not produce biological findings.
