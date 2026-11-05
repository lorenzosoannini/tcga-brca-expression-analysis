# TCGA Breast Cancer Gene Expression Analysis

**From tumor–normal expression differences to the biological interpretation of breast cancer pathways.**

This Bioinformatics exam project investigates how gene expression differs between
breast invasive carcinoma (BRCA) and matched normal tissue, which biological
processes are associated with those differences, and how selected genes help
interpret the findings.

The analysis moves from a broad expression matrix to a selected gene set, explores
its structure with clustering and PCA, then uses functional enrichment and
literature to connect the statistical results to biological questions.

**R · TCGA-BRCA · Differential expression · PCA · Functional enrichment**

**[Read the original exam presentation](docs/Bioinformatics%20-%20BRCA.pdf)** —
the reference for the study's narrative, figures and biological discussion.
All exam findings discussed below are already available in the presentation and
bundled result exports. Reading or reviewing this project requires no R execution.

## The analysis, step by step

### 1. Build a paired tumor–normal comparison

*Presentation: slides 4–7.*

The starting TCGA-BRCA matrix contains **56,602 genes and 1,222 samples**.
Rather than comparing unrelated groups, the analysis retains the **112 patients
with both tumor and normal samples**, giving 224 samples in a paired design.

Genes with zero overall expression are removed, expression is transformed with
`log2(x + 1)`, and low-variability genes are excluded using an IQR filter.
This leaves **32,456 genes** for the next selection stage.

### 2. Identify the strongest expression differences

*Presentation: slides 8–13.*

The next step combines an expression-difference threshold with paired t-tests
and adjusted p-values. The resulting set contains **1,556 differentially expressed
genes: 513 upregulated and 1,043 downregulated** in tumor tissue.

Volcano plots summarize the selected differences; boxplots highlight **COL10A1**,
the most upregulated gene, and **FABP4**, the most downregulated. The heatmap then
asks whether the selected expression patterns organize the samples by condition.

![Original heatmap of the 1,556 selected genes](docs/assets/heatmap.png)

*The selected genes show contrasting expression patterns between tumor and
normal tissue. Rows are standardized; clustering uses correlation distance
and average linkage.*

Outputs: [gene statistics](results/exam/BRCA/DEG.txt),
[volcano plot](results/exam/BRCA/volcano_plot_pval.pdf),
[COL10A1 boxplot](results/exam/BRCA/boxplot_max.pdf),
[FABP4 boxplot](results/exam/BRCA/boxplot_min.pdf).

### 3. Explore the structure of the selected expression profiles

*Presentation: slides 14–17.*

PCA provides a complementary view of the same selected gene set: how much of its
variation can be summarized in a few dimensions, and which genes contribute to
those dimensions?

**PC1 explains 50.9% and PC2 10.1% of the variance**, for 61.0% together. The score
plot shows the tumor and normal groups separating in this representation. The
contribution analysis highlights **LINC02202 and CAV1**, which are revisited in
the biological discussion.

![Original PCA score plot](docs/assets/pca.png)

*PCA describes genes already selected using the tumor/normal labels. This
separation is exploratory and is not an independent classification evaluation.*

Outputs: [variance plots](results/exam/PCA/pareto_scree_plots.pdf),
[gene contributions](results/exam/PCA/contributions_plot.pdf).

### 4. Move from individual genes to biological processes

*Presentation: slides 18–22.*

Differential expression identifies genes, but interpreting a list of 1,556 genes
requires a broader view. UP and DOWN lists are therefore analyzed separately
with Enrichr, using GO, KEGG and DisGeNET annotations.

The presentation focuses on two complementary views:

| Analysis | Main findings | Role in the study |
| --- | --- | --- |
| KEGG, upregulated genes | **Cell cycle** is the leading term by adjusted p-value | Connect the UP list to processes discussed in relation to proliferation |
| KEGG, downregulated genes | **Regulation of lipolysis in adipocytes**, **PPAR signaling pathway**, **Axon guidance** | Guide the discussion toward lipid metabolism, tissue composition and signaling |
| DisGeNET, both directions | Breast carcinoma, mammary neoplasms and metastasis-related annotations | Place the gene lists in a disease-annotation context |

![Original KEGG enrichment plot for upregulated genes](docs/assets/kegg-up.png)

*Archived KEGG 2021 Human output. GO outputs are also retained in the repository;
the slides emphasize KEGG and DisGeNET.*

Outputs: [UP enrichment tables and figures](results/exam/BRCA/Functional_Enrichment_UP),
[DOWN enrichment tables and figures](results/exam/BRCA/Functional_Enrichment_DOWN).

### 5. Return to selected genes and interpret the findings

*Presentation: slides 23–32.*

The final part brings together three sources of evidence: genes with the largest
expression changes, genes contributing to PCA, and genes connected to enriched
pathways. The literature discussion then examines what those observations might
mean and where a simple interpretation would be insufficient.

| Gene or theme | Why it is discussed | Interpretation explored in the presentation |
| --- | --- | --- |
| **COL10A1** | Largest positive expression difference | Relate the observed overexpression to published breast cancer findings |
| **FABP4** | Largest negative expression difference | Examine the apparent tension with literature on adipocyte–tumor interactions; consider differences in adipocyte content as a possible explanation |
| **LINC02202 and CAV1** | Highlighted by PCA contribution analysis | Discuss their biological context and, for CAV1, its context-dependent roles |
| **CDK1, cell cycle and p53 signaling** | Connect the gene-level results to the pathway discussion | Examine the relationship between proliferation and cell-cycle regulation; the dataset does not provide statistical evidence of differential TP53 expression under the chosen procedure |

The FABP4 example is particularly useful: a strong statistical difference does not
by itself explain the underlying biology. Tissue composition is discussed as a
hypothesis, rather than treated as a mechanism established by these data.
Likewise, unchanged TP53 expression under this analysis does not establish normal
p53 pathway function.

## Scope and methodological notes

This is an **exploratory analysis of a paired cohort**, with interpretation informed
by pathway annotations and literature. It does not include an independently
validated diagnostic model or a BRCA survival study.

The original selection order is preserved: IQR filtering, an absolute difference
of at least 1 between mean log2-transformed expression values, then paired tests
and Benjamini–Hochberg adjustment at 0.01. Fold-change filtering before adjustment
affects the multiple-testing interpretation; PCA uses the already-selected
genes. Input processing metadata were not fully recorded.
See [methodology and limitations](docs/METHODOLOGY.md) for the details.

## Explore the repository

| Location | Contents |
| --- | --- |
| [Original presentation](docs/Bioinformatics%20-%20BRCA.pdf) | Full exam narrative, plots, interpretation and bibliography |
| [BRCA results](results/exam/BRCA) | Gene statistics, selected expression matrix, sample lists, enrichment outputs and figures |
| [PCA results](results/exam/PCA) | Original score, variance and contribution plots |
| [Original R scripts](archive/code) | Unmodified core analysis scripts from the exam |
| [Portable R scripts](scripts) | Optional versions with portable paths and input checks |
| [Provenance](docs/PROVENANCE.md) | Artifact origins and the changes made for this portfolio copy |

The original presentation, tables and PDF figures are preserved byte-for-byte.
Lightweight automated checks verify their integrity and internal consistency;
they do not rerun the R analyses. See [validation scope](docs/VALIDATION.md).

For readers who want to run the code themselves, [optional execution instructions](docs/RUNNING.md)
cover dependencies and inputs. R execution is not needed to inspect the results.
The large source matrix and clinical files are excluded from GitHub.

## References and attribution

The original presentation contains the bibliography used in the biological
discussion. Data are TCGA-derived; enrichment annotations come from the named
databases. Course scaffolding and third-party material retain their provenance,
as described in [the attribution notes](docs/PROVENANCE.md).

Project by **Lorenzo Soannini**, originally developed for the Bioinformatics exam
in 2024 and organized for this portfolio in 2026.
