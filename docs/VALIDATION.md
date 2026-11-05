# Validation scope

The project presents existing exam results. Its automated workflow checks those
artifacts without executing R or recomputing the statistical analyses.

`python3 scripts/verify_results.py` checks SHA-256 hashes of copied exam artifacts;
the 1,556-gene total and UP/DOWN counts; finite values and selection thresholds;
patient pairing and matrix order; consistency of exported mean log-expression
differences with the saved matrix; extreme genes; and leading KEGG terms.
It uses only the Python standard library and does not query external databases.

During initial preparation, a one-time R 4.4.3 check successfully executed the
portable PCA on the bundled matrix, verified PC1 and PC2 against the slides, and
exercised the portable expression/PCA scripts on synthetic paired data. This
checked the optional scripts, not the biological findings. No new full-cohort
differential expression analysis or live enrichment query was performed.

The original exam presentation and outputs remain the reference for reported
results. The historical package versions and exact source-data release were not
recovered. Optional execution instructions are in RUNNING.md.
