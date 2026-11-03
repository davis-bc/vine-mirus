"""
mock_cenote_taker3.py — synthetic Cenote-Taker3 virus summary for test_mode.
Called via shell: python mock_cenote_taker3.py <contigs> <output> <min_len>

Reproduces the real quirks the evidence parser must handle:
  - "contig" is Cenote-Taker3's OWN internal ID (cenote_taker3_N) and cannot be
    mapped back to the input; "input_name" is the one to key on
  - "input_name" holds the FULL original FASTA header, not the bare ID
"""

import argparse
import random
from pathlib import Path

import pandas as pd
from Bio import SeqIO

random.seed(2909)

parser = argparse.ArgumentParser()
parser.add_argument("contigs")
parser.add_argument("output")
parser.add_argument("min_len", type=int)
args = parser.parse_args()

LINEAGES = [
    "Viruses;Riboviria;Orthornavirae;Pisuviricota;Pisoniviricetes;Picornavirales;Marnaviridae",
    "Viruses;Monodnaviria;Shotokuvirae;Cressdnaviricota;Repensiviricetes;Geplafuvirales;Genomoviridae",
    "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes",
    "Viruses;Riboviria;Orthornavirae;Kitrinoviricota;Tolucaviricetes;Tolivirales;Tombusviridae",
]

contigs = [r for r in SeqIO.parse(args.contigs, "fasta") if len(r.seq) >= args.min_len]
Path(args.output).parent.mkdir(parents=True, exist_ok=True)

rows = []
for i, rec in enumerate(contigs):
    if random.random() > 0.25:      # reports only annotated virus hits
        continue
    rows.append({
        "contig": f"cenote_taker3_{i}",   # internal ID: NOT usable for joining
        "input_name": f"{rec.id} flag=0 multi={random.uniform(2, 9):.4f} len={len(rec.seq)}",
        "virus_seq_length": len(rec.seq),
        "end_feature": random.choice(["DTR", "ITR", "None"]),
        "taxonomy_hierarchy": random.choice(LINEAGES),
        "organism_name": "virus sp.",
    })

cols = ["contig", "input_name", "virus_seq_length", "end_feature", "taxonomy_hierarchy", "organism_name"]
pd.DataFrame(rows, columns=cols).to_csv(args.output, sep="\t", index=False)
print(f"[mock] Cenote-Taker3: {len(rows)}/{len(contigs)} contigs annotated")
