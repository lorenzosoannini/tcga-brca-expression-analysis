# Methodology and interpretation

## Historical analysis

The reference for the exam findings is [the original presentation](Bioinformatics%20-%20BRCA.pdf).
The available R scripts and result exports accompany that account. The portable
scripts retain the selection procedure; they do not replace it with a new study.

1. Select normal (`1x`) and tumor (`0x`) TCGA sample types, keep the first sample
   per patient and condition, and align both groups by patient identifier.
2. Retain the 112 patients with both conditions; exclude genes with overall zero
   mean and transform expression with `log2(x + 1)`.
3. Exclude genes at or below the 40th percentile of IQR across the paired cohort.
   The presentation reports 54,094 genes after zero filtering and 32,456 after IQR.
4. Compute the difference between the tumor and normal **means of log-transformed
   values**. Retain absolute differences of at least 1 (`thr_fc = 2` in the original
   script). This is the code's definition; it is not the log ratio of arithmetic
   means written on slide 8.
5. Apply a paired t-test to each remaining gene and Benjamini-Hochberg adjustment
   within this prefiltered set. Retain adjusted p-values at or below 0.01.
   The slides report 1,557 genes before this final filter and 1,556 after it.
6. Cluster the selected expression matrix using correlation distance and average
   linkage, with rows standardized for the heatmap. Apply centered, scaled PCA
   to this same selected matrix.
7. Query UP and DOWN lists separately with Enrichr against GO Molecular Function
   2021, GO Biological Process 2021, KEGG 2021 Human and DisGeNET. Archived
   enrichment tables retain terms with adjusted p-values below 0.05.

## Limits that affect interpretation

- **Selection before multiple testing:** fold-change filtering uses the same
  condition contrast as the tests. The adjusted p-values belong to this historical
  prefiltered test family; they should not be advertised as a demonstrated
  genome-wide FDR guarantee. A future analysis should reconsider the testing
  procedure after establishing the input units and processing provenance.
- **PCA after supervised feature selection:** the genes were selected using the
  tumor/normal labels before PCA. The visible separation is a description of the
  selected cohort, not independent evidence of predictive performance. There is
  no held-out classifier, cross-validation or external cohort in this project.
- **Paired tests:** 112 pairs do not automatically guarantee normality. The
  original script's sample-size comment is an oversimplification. Diagnostics
  should concern the paired differences and any influential observations.
- **Biological interpretation:** tissue composition and other confounders can
  affect bulk expression. The FABP4 discussion in the slides is a hypothesis,
  not a mechanism established by this analysis. Enrichment associations do not
  establish causality, infection status or therapeutic efficacy.
- **Enrichment background:** the historical Enrichr call did not specify the
  measured/tested gene universe. This can affect over-representation results.
  Live database contents and term membership can change over time.
- **Scope:** this is an educational exploratory analysis. The survival scripts in
  the source folder use example lung/myeloma data and do not establish BRCA
  survival associations. They are excluded from this repository.

## Portable scripts

The cleaned expression script uses the same thresholds and pairing convention,
adds input checks, fixes empty-index subsetting, and writes to a separate output
directory. Its plots are updated visualizations of a rerun. The PCA script uses
`prcomp(center = TRUE, scale. = TRUE)` as before, with ggplot2 instead of the
original ggbiplot-based plotting code. Enrichment saves full live responses and
significant subsets for both directions. Every script records `sessionInfo()`.

The complete historical environment was not saved. Dependencies are installed
from CRAN; the CI environment is documented in its workflow, not represented as
an exact recreation of the 2024 package versions. Historical figures and tables
are never overwritten by these commands.
