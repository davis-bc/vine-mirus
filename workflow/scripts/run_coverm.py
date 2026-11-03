"""
run_coverm.py — map the QC'd reads of every sample against the abundance
reference (vOTU representatives when clustering is on, otherwise the viral
contigs) and write one abundance matrix per metric.

Mapping uses minimap2-sr: short-read mapping against a contig set this size is
minutes with minimap2 versus roughly an hour of BWA indexing, and CoverM drives
it internally so no separate alignment rule is needed.

Reads are the per-sample QC'd (optionally host-removed/deduplicated) reads, NOT
the bbnorm-normalized ones — normalization distorts depth by design and would
make abundance meaningless.

Two outputs:
  - abundance_matrix.tsv      the primary metric (config abundance.primary_metric),
                              one column per sample, bare sample names
  - abundance_all_metrics.tsv every requested metric, CoverM's own column names

Called as a Snakemake script.
"""

import random
import sys
from pathlib import Path

import pandas as pd

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

params = snakemake.params
test_mode = bool(params.test_mode)
samples = list(params.samples)
primary = str(params.primary_metric)
metrics = list(params.metrics)          # primary first, then the extras

primary_out = snakemake.output.abundance
all_out = snakemake.output.all_metrics
Path(primary_out).parent.mkdir(parents=True, exist_ok=True)

# CoverM's column label for each method it can report.
METRIC_LABEL = {
    "mean": "Mean",
    "trimmed_mean": "Trimmed Mean",
    "covered_bases": "Covered Bases",
    "covered_fraction": "Covered Fraction",
    "variance": "Variance",
    "length": "Length",
    "count": "Read Count",
    "reads_per_base": "Reads per base",
    "rpkm": "RPKM",
    "tpm": "TPM",
    "relative_abundance": "Relative Abundance (%)",
}

reference = snakemake.input.reference
contig_ids = [line[1:].strip().split()[0]
              for line in open(reference) if line.startswith(">")]
print(f"Abundance reference: {reference} ({len(contig_ids)} sequences)")
print(f"Samples: {len(samples)}; metrics: {metrics} (primary: {primary})")

if not contig_ids:
    raise SystemExit(f"Abundance reference {reference} contains no sequences")

# ── Run CoverM (or mock it) ───────────────────────────────────────────────────
if test_mode:
    print("[mock] CoverM")
    random.seed(42)
    rows = []
    for contig in contig_ids:
        row = {"Contig": contig}
        for sample in samples:
            for metric in metrics:
                label = METRIC_LABEL.get(metric, metric.upper())
                row[f"{sample}_R1.fastq.gz {label}"] = round(random.uniform(0.0, 1000.0), 2)
        rows.append(row)
    raw = pd.DataFrame(rows)
else:
    from snakemake.shell import shell

    coupled = " ".join(f"{r1} {r2}" for r1, r2 in zip(snakemake.input.r1,
                                                      snakemake.input.r2))
    # One --methods flag taking every metric, so the mapping happens once.
    methods = "--methods " + " ".join(metrics)
    shell(
        "coverm contig"
        " --reference {reference}"
        f" --coupled {coupled}"
        " --mapper minimap2-sr"
        f" {methods}"
        " --threads {snakemake.threads}"
        " --output-file {all_out}"
        " >> {snakemake.log[0]} 2>&1"
    )
    raw = pd.read_csv(all_out, sep="\t")

if test_mode:
    raw.to_csv(all_out, sep="\t", index=False)

print(f"\nCoverM table: {raw.shape[0]} rows x {raw.shape[1]} columns")
print("Columns: " + ", ".join(map(str, raw.columns)))


# ── Reduce the primary metric to one bare-sample-name column each ─────────────
def sample_name(column):
    """Strip CoverM's read-set decoration down to the sample name."""
    name = str(column).split()[0].split("/")[-1]
    for suffix in ("_R1.fastq.gz", "_R2.fastq.gz", "_R1.fq.gz", "_R1", ".fastq.gz",
                   ".fq.gz", ".fa", ".fasta"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


primary_label = METRIC_LABEL.get(primary, primary.upper()).lower()
# Every label this run could have produced, longest first: "Mean" is a suffix
# of "Trimmed Mean", so a plain endswith() would claim both columns for `mean`
# and one sample's value would silently overwrite the other.
all_labels = sorted((METRIC_LABEL.get(m, m.upper()).lower() for m in metrics),
                    key=len, reverse=True)

contig_col = raw.columns[0]
matrix = {"contig_id": raw[contig_col]}
for column in raw.columns[1:]:
    text = str(column)
    match = next((lab for lab in all_labels if text.lower().endswith(lab)), None)
    if match != primary_label:
        continue
    # Everything before the label is CoverM's read-set name.
    stem = text[: len(text) - len(match)].strip() or text
    matrix[sample_name(stem)] = raw[column]

df = pd.DataFrame(matrix)
found = [c for c in df.columns[1:]]
print(f"\nPrimary metric ({primary}) columns: {found}")

missing = [s for s in samples if s not in found]
unexpected = [c for c in found if c not in samples]
if missing or unexpected:
    # A silently mis-named column would produce an all-NaN abundance column in
    # the summary, which is exactly the kind of clean-exit failure to catch here.
    raise SystemExit(
        f"CoverM columns do not match the sample list. Missing: {missing}; "
        f"unexpected: {unexpected}"
    )

df.to_csv(primary_out, sep="\t", index=False)
print(f"\nPrimary matrix  : {df.shape[0]} contigs x {len(found)} samples -> {primary_out}")
print(f"All metrics     : {raw.shape[1] - 1} columns -> {all_out}")

totals = df[found].fillna(0)
print(f"Column sums     : "
      + ", ".join(f"{c}={totals[c].sum():.2f}" for c in found))
if totals.to_numpy().sum() == 0:
    print("[WARN] every abundance value is zero: check that the reads actually map")
log.close()
