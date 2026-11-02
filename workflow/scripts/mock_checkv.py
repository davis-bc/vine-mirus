"""
mock_checkv.py — synthetic CheckV quality_summary.tsv for test_mode.
Called via shell: python mock_checkv.py <candidates_fasta> <output>

Column names match the real quality_summary.tsv, including contig_id.
"""

import argparse
import random
from pathlib import Path

import pandas as pd
from Bio import SeqIO

random.seed(7331)

parser = argparse.ArgumentParser()
parser.add_argument("candidates")
parser.add_argument("output")
args = parser.parse_args()

QUALITIES = ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]

rows = []
for rec in SeqIO.parse(args.candidates, "fasta"):
    completeness = round(random.uniform(5, 100), 2)
    rows.append({
        "contig_id": rec.id,
        "contig_length": len(rec.seq),
        "provirus": "No",
        "proviral_length": "NA",
        "gene_count": random.randint(1, 25),
        "viral_genes": random.randint(0, 12),
        "host_genes": random.randint(0, 3),
        "checkv_quality": random.choice(QUALITIES),
        "miuvig_quality": "Genome-fragment",
        "completeness": completeness,
        "completeness_method": "AAI-based",
        "contamination": round(random.uniform(0, 5), 2),
        "kmer_freq": 1.0,
        "warnings": "",
    })

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
print(f"[mock] CheckV: {len(rows)} contigs -> {args.output}")
