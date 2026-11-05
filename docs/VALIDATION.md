# Validation scope

`python3 scripts/verify_results.py` checks the SHA-256 hashes of copied exam
artifacts; the 1,556-gene total and UP/DOWN counts; finite values and selection
thresholds; patient pairing and matrix order; every exported mean log-expression
difference against the selected matrix; extreme genes; and leading KEGG terms.
This check uses only the Python standard library.

The GitHub Actions workflow additionally parses the portable R files, runs PCA
on the original selected matrix, checks PC1 and PC2 against the presentation's
rounded values, and executes expression analysis plus PCA on synthetic paired
data with shuffled input columns and known up/down effects. Synthetic output is
used solely to test execution and matching, never as a biological result.

Full differential expression on the large historical input and live enrichment
queries are not run in CI. The former requires the locally retained matrix; the
latter depends on external database/service state. The historical package
versions and exact source-data release were not recovered.
