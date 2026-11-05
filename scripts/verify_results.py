"""Validate archived exam results without requiring R or the full source matrix."""
import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAM = ROOT / "results/exam/BRCA"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    hashes = json.loads((ROOT / "docs/artifact-checksums.json").read_text())
    for name, expected in hashes.items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        require(actual == expected, f"Historical artifact changed: {name}")
    with (EXAM / "DEG.txt").open() as handle:
        genes = list(csv.DictReader(handle, delimiter="\t"))
    require(len(genes) == 1556, "Expected 1,556 selected genes")
    require(Counter(g["direction"] for g in genes) == {"UP": 513, "DOWN": 1043}, "Direction counts differ")
    require(len({g["ensemblID"] for g in genes}) == 1556, "Duplicate gene IDs")
    for g in genes:
        fc, p, q = map(float, (g["logFC"], g["pvalue"], g["pval_adj"]))
        require(all(map(math.isfinite, (fc, p, q))), "Non-finite gene statistics")
        require(abs(fc) >= 1 and 0 <= p <= q <= 0.01, f"Selection criteria differ: {g['genes']}")
        require(g["direction"] == ("UP" if fc > 0 else "DOWN"), "Direction/sign mismatch")
    require(max(genes, key=lambda g: float(g["logFC"]))["genes"] == "COL10A1", "UP extreme differs")
    require(min(genes, key=lambda g: float(g["logFC"]))["genes"] == "FABP4", "DOWN extreme differs")
    normal = (EXAM / "list_normal.txt").read_text().splitlines()
    tumor = (EXAM / "list_cancer.txt").read_text().splitlines()
    require(len(normal) == len(tumor) == 112, "Expected 112 pairs")
    def patient(sample):
        return "-".join(sample.split("-")[:3])
    require(len(set(map(patient, normal))) == 112, "Repeated patient")
    require(list(map(patient, normal)) == list(map(patient, tumor)), "Pair ordering differs")
    for samples, kind in [(normal, "1"), (tumor, "0")]:
        require(all(re.fullmatch(r"TCGA-\w+-\w+-" + kind + r"\d.*", s) for s in samples), "Sample type differs")
    with (EXAM / "matrix_DEG.txt").open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        require(header[1:] == normal + tumor, "Matrix sample order differs")
        matrix = list(reader)
    require(len(matrix) == len(genes), "Gene matrix dimensions differ")
    for row, gene in zip(matrix, genes):
        require(row[0] == gene["genes"] + "|" + gene["ensemblID"], "Gene order differs")
        require(len(row) == 225, "Expected 224 expression values")
        values = list(map(float, row[1:]))
        require(all(math.isfinite(v) and v >= 0 for v in values), "Invalid log expression")
        fc = sum(values[112:]) / 112 - sum(values[:112]) / 112
        require(math.isclose(fc, float(gene["logFC"]), abs_tol=1e-10), "Expression/FC mismatch")
    for direction, expected in [("UP", "Cell cycle"), ("DOWN", "Regulation of lipolysis in adipocytes")]:
        with (EXAM / f"Functional_Enrichment_{direction}/KEGG_adj_pval_0.05.txt").open() as handle:
            terms = list(csv.DictReader(handle, delimiter="\t"))
        require(min(terms, key=lambda r: float(r["Adjusted.P.value"]))["Term"] == expected, "Leading KEGG term differs")
    print(f"PASS: {len(hashes)} unchanged historical artifacts; 112 aligned pairs; 1,556 genes; expression differences; KEGG leaders.")


if __name__ == "__main__":
    main()
