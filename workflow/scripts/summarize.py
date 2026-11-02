"""
summarize.py — join the evidence table with the CoverM abundance matrix into
the final per-contig summary: one row per viral contig, one abundance column
per sample.

Only contigs present in viral_contigs.fa appear here (candidates, with
high-confidence phage excluded unless phage_filtering.include_phage). The full
evidence table, including everything that was screened out, stays in
results/evidence_table.tsv.

Called as a Snakemake script.
"""

import sys
from pathlib import Path

import pandas as pd

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

samples = list(snakemake.params.samples)
metric = snakemake.params.metric

# ── Contigs that made it into the viral FASTA ─────────────────────────────────
viral_ids = []
with open(snakemake.input.viral_fa) as fh:
    for line in fh:
        if line.startswith(">"):
            viral_ids.append(line[1:].strip().split()[0])
print(f"Viral contigs: {len(viral_ids)}")

evidence = pd.read_csv(snakemake.input.evidence, sep="\t")
df = evidence[evidence["contig_id"].isin(viral_ids)].copy()
if len(df) != len(viral_ids):
    missing = sorted(set(viral_ids) - set(df["contig_id"]))[:10]
    raise SystemExit(f"{len(viral_ids) - len(df)} FASTA contigs are absent from the "
                     f"evidence table, e.g. {missing}")

# ── Abundance ─────────────────────────────────────────────────────────────────
abundance = pd.read_csv(snakemake.input.abundance, sep="\t")
# CoverM writes "Contig" then one column per BAM/read set, named after the read
# file with the metric appended ("sample1_R1.fastq.gz RPKM"). Reduce each to the
# bare sample name so the columns are predictable.
renamed = ["contig_id"]
for col in abundance.columns[1:]:
    name = col.split()[0].split("/")[-1]
    for suffix in ("_R1.fastq.gz", "_R2.fastq.gz", "_R1", ".fastq.gz", ".fa", ".fasta"):
        name = name.replace(suffix, "")
    renamed.append(name)
abundance.columns = renamed

unexpected = [c for c in abundance.columns[1:] if c not in samples]
if unexpected:
    print(f"[WARN] abundance columns not matching a sample name: {unexpected}")
missing_samples = [s for s in samples if s not in abundance.columns]
if missing_samples:
    print(f"[WARN] no abundance column for: {missing_samples}")

df = df.merge(abundance, on="contig_id", how="left")

dups = df["contig_id"].duplicated(keep=False)
if dups.any():
    raise SystemExit("Duplicate contig_id rows after merge: "
                     f"{df.loc[dups, 'contig_id'].unique()[:10].tolist()}")

# ── Column order ──────────────────────────────────────────────────────────────
meta = [c for c in ["contig_id", "sample_of_origin", "length", "viral_score",
                    "n_tools_called", "confidence_tier", "taxonomy", "host_domain",
                    "checkv_quality", "checkv_completeness"] if c in df.columns]
abundance_cols = [c for c in df.columns if c in samples]
df = df[meta + abundance_cols].sort_values(["confidence_tier", "viral_score"],
                                           ascending=[True, False])

out = snakemake.output.tsv
Path(out).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, sep="\t", index=False)

print(f"Summary: {len(df)} contigs x {len(abundance_cols)} sample {metric.upper()} columns -> {out}")
if abundance_cols:
    print("Abundance column sums:")
    for col in abundance_cols:
        print(f"  {col}: {df[col].sum():.2f}")
    if df[abundance_cols].fillna(0).to_numpy().sum() == 0:
        print("[WARN] every abundance value is zero or missing - check read mapping")
log.close()
