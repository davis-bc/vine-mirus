"""
summarize.py — join the evidence table, vOTU membership and the CoverM
abundance matrix into the final per-contig summary.

Only contigs present in viral_contigs.fa appear here (candidates, with
high-confidence phage excluded unless phage_filtering.include_phage). The full
evidence table, including everything that was screened out, stays in
results/evidence_table.tsv.

Abundance is measured per vOTU: the reads of a cluster are mapped to its
representative, so every member row reports its vOTU's abundance and
is_votu_representative marks which contig carried the mapping. With clustering
disabled, abundance is per contig.

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

# ── vOTU membership ───────────────────────────────────────────────────────────
# votu_cluster.py clusters viral_contigs.fa itself, so every row here must be
# present; a gap means the two steps saw different contig sets.
clusters = getattr(snakemake.input, "clusters", None)
if clusters:
    votu = pd.read_csv(clusters, sep="\t")
    df = df.merge(votu, on="contig_id", how="left")
    orphans = df.loc[df["votu_representative"].isna(), "contig_id"].tolist()
    if orphans:
        raise SystemExit(f"{len(orphans)} viral contigs are missing from the vOTU "
                         f"cluster table, e.g. {orphans[:10]}")
    n_reps = int(df["is_votu_representative"].sum())
    print(f"vOTUs: {n_reps} representatives for {len(df)} contigs")
    abundance_key = "votu_representative"
else:
    print("vOTU clustering disabled: abundance is per contig")
    abundance_key = "contig_id"

# ── Abundance ─────────────────────────────────────────────────────────────────
abundance = pd.read_csv(snakemake.input.abundance, sep="\t")
# run_coverm.py already reduced the columns to bare sample names; verify rather
# than assume, since a mismatch here silently empties the abundance columns.
if abundance.columns[0] != "contig_id":
    raise SystemExit(f"abundance matrix's first column is {abundance.columns[0]!r}, "
                     "expected 'contig_id'")
found = list(abundance.columns[1:])
missing_samples = [s for s in samples if s not in found]
if missing_samples:
    print(f"[WARN] no abundance column for: {missing_samples}")
unexpected = [c for c in found if c not in samples]
if unexpected:
    print(f"[WARN] abundance columns not matching a sample name: {unexpected}")

if found:
    mapped = set(abundance["contig_id"])
    keys = set(df[abundance_key].dropna())
    unmeasured = sorted(keys - mapped)
    if unmeasured:
        print(f"[WARN] {len(unmeasured)} of {len(keys)} abundance references have no "
              f"row in the matrix, e.g. {unmeasured[:5]}")

df = df.merge(abundance.rename(columns={"contig_id": abundance_key}),
              on=abundance_key, how="left")

dups = df["contig_id"].duplicated(keep=False)
if dups.any():
    raise SystemExit("Duplicate contig_id rows after merge: "
                     f"{df.loc[dups, 'contig_id'].unique()[:10].tolist()}")

# ── Column order ──────────────────────────────────────────────────────────────
meta = [c for c in ["contig_id", "sample_of_origin", "length", "viral_score",
                    "n_tools_called", "confidence_tier", "taxonomy", "host_domain",
                    "checkv_quality", "checkv_completeness",
                    "votu_representative", "votu_cluster_size",
                    "is_votu_representative"] if c in df.columns]
abundance_cols = [c for c in df.columns if c in samples]
sort_by = ["confidence_tier", "viral_score"]
df = df[meta + abundance_cols].sort_values(sort_by, ascending=[True, False])

out = snakemake.output.tsv
Path(out).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, sep="\t", index=False)

print(f"\nSummary: {len(df)} contigs x {len(abundance_cols)} sample "
      f"{metric.upper()} columns -> {out}")
if abundance_cols:
    note = (" (per-contig rows, so vOTU members repeat their representative)"
            if clusters else "")
    print(f"Abundance column sums{note}:")
    for col in abundance_cols:
        print(f"  {col}: {df[col].sum():.2f}")
    if df[abundance_cols].fillna(0).to_numpy().sum() == 0:
        print("[WARN] every abundance value is zero or missing - check read mapping")
log.close()
