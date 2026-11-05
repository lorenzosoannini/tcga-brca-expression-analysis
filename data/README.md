# Input data

The exam used a TCGA-BRCA expression matrix with 56,602 gene rows and 1,222 sample
columns (113 normal and 1,109 tumor samples), yielding 112 matched patient pairs.
These dimensions are reported in the original presentation, slide 4.

The 635,695,835-byte source file `matrix_RNAseq_BRCA.txt` and the separate clinical
and intermediate CSV files are excluded from Git. The original files remain in
the local exam folder. To run the expression workflow, copy the original matrix
here, or pass its path as the first command-line argument.

Expected format: tab-separated, header row with TCGA sample barcodes, first column
with unique `geneSymbol|ensemblID` identifiers, remaining columns with finite,
non-negative expression values on the original **unlogged** input scale.
The script performs `log2(x + 1)` itself. Do not supply raw read counts or an
already-log-transformed replacement and assume it reproduces the exam results.

The exact TCGA processing release, normalization unit and original download URL
were not recorded in the available project files. For that reason this repository
does not provide a fabricated download command or claim that a current TCGA
matrix reproduces the historical input. TCGA data access information is available
from the [NCI Genomic Data Commons](https://gdc.cancer.gov/).

The derived selected-gene matrix and sample lists under `results/exam/BRCA/` are
included so PCA can run without the large input file. These are historical
TCGA-derived exports, not synthetic examples. No clinical table is bundled.
