"""
summarize.py — join viral scores, taxonomy, CheckV quality, and CoverM abundance
into the final viral_contigs_abundance.tsv. Called as a Snakemake script.

Output columns:
  contig_id, sample_of_origin, length, viral_score, taxonomy,
  checkv_quality, <sample1_rpkm>, <sample2_rpkm>, ...
"""

import sys
from pathlib import Path

import pandas as pd

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

scores_path    = snakemake.input.scores
taxonomy_path  = snakemake.input.taxonomy
checkv_path    = snakemake.input.checkv
abundance_path = snakemake.input.abundance
out_path       = snakemake.output.tsv

# ── Load inputs ───────────────────────────────────────────────────────────────
scores   = pd.read_csv(scores_path, sep="\t")
passing  = scores[scores["PASS"]][["sample", "contig_id", "viral_score"]]

taxonomy = pd.read_csv(taxonomy_path, sep="\t")[["contig_id", "sample", "taxonomy"]]

checkv   = pd.read_csv(checkv_path, sep="\t")
checkv   = checkv[["contig_name", "contig_length", "checkv_quality"]].rename(
    columns={"contig_name": "contig_id", "contig_length": "length"}
)

abundance = pd.read_csv(abundance_path, sep="\t")
# CoverM output: first column = Contig, then one column per sample
# Normalise column names: strip coverm path/suffix noise
abundance.columns = [
    "contig_id" if i == 0 else c.split(" ")[0].split("/")[-1].replace("_R1.fastq.gz", "")
    for i, c in enumerate(abundance.columns)
]

# ── Merge ─────────────────────────────────────────────────────────────────────
df = passing.merge(taxonomy, on=["contig_id", "sample"], how="left")
df = df.merge(checkv,       on="contig_id",              how="left")
df = df.merge(abundance,    on="contig_id",              how="left")

df = df.rename(columns={"sample": "sample_of_origin"})

# Reorder: metadata first, then per-sample RPKM columns
meta_cols   = ["contig_id", "sample_of_origin", "length", "viral_score", "taxonomy", "checkv_quality"]
sample_cols = [c for c in df.columns if c not in meta_cols]
df = df[meta_cols + sample_cols]

# contig_id is globally unique ("{sample}__<id>", see prep_contigs.py), so a
# duplicate here means a join fanned out. Fail instead of silently dropping rows.
dups = df["contig_id"].duplicated(keep=False)
if dups.any():
    raise ValueError(f"Duplicate contig_id rows after merge: {df.loc[dups, 'contig_id'].unique()[:10].tolist()}")

Path(out_path).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out_path, sep="\t", index=False)
print(f"Final summary: {len(df)} viral contigs written to {out_path}")
log.close()
