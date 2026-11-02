"""
mock_vibrant.py — synthetic VIBRANT genome-quality table for test_mode.
Called via shell: python mock_vibrant.py <contigs> <output> <min_len>

Deliberately reproduces the two real quirks the evidence parser must handle:
  - "scaffold" holds the FULL original FASTA header, not the bare contig ID
  - the quality column is "Quality" (capital Q)
VIBRANT lists only contigs it calls viral, so a subset is emitted.
"""

import argparse
import random
from pathlib import Path

import pandas as pd
from Bio import SeqIO

random.seed(1301)   # distinct per tool, so tools disagree as they would in reality

parser = argparse.ArgumentParser()
parser.add_argument("contigs")
parser.add_argument("output")
parser.add_argument("min_len", type=int)
args = parser.parse_args()

QUALITIES = ["complete circular", "high quality draft", "medium quality draft", "low quality draft"]

contigs = [r for r in SeqIO.parse(args.contigs, "fasta") if len(r.seq) >= args.min_len]
Path(args.output).parent.mkdir(parents=True, exist_ok=True)

rows = []
for rec in contigs:
    if random.random() > 0.35:      # VIBRANT reports only its viral calls
        continue
    rows.append({
        # full header, exactly as VIBRANT writes it
        "scaffold": f"{rec.id} flag=1 multi={random.uniform(2, 9):.4f} len={len(rec.seq)}",
        "type": random.choice(["lytic", "lysogenic"]),
        "Quality": random.choice(QUALITIES),
    })

pd.DataFrame(rows, columns=["scaffold", "type", "Quality"]).to_csv(args.output, sep="\t", index=False)
print(f"[mock] VIBRANT: {len(rows)}/{len(contigs)} contigs called")
