"""
split_fasta.py — split a contig FASTA into N chunks of roughly equal total
bases, so the slow per-contig tools (geNomad, VirSorter2, VIBRANT,
Cenote-Taker3) can run as parallel jobs. Called as a Snakemake script.

Balancing by bases, not by contig count, matters because contig length
distributions are heavily skewed: a chunk of few long contigs otherwise
dominates the wall time.

Every requested chunk file is always written, empty if there is nothing left
to put in it, so the DAG never has a missing output.
"""

import sys
from pathlib import Path

from Bio import SeqIO

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

outputs = list(snakemake.output.chunks)
n = len(outputs)

records = sorted(SeqIO.parse(snakemake.input.contigs, "fasta"),
                 key=lambda r: len(r.seq), reverse=True)

# Greedy longest-first bin packing: each contig goes to the lightest chunk.
bins = [[] for _ in range(n)]
loads = [0] * n
for rec in records:
    i = loads.index(min(loads))
    bins[i].append(rec)
    loads[i] += len(rec.seq)

for path, bucket in zip(outputs, bins):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for rec in bucket:
            fh.write(f">{rec.id}\n{str(rec.seq)}\n")

print(f"Split {len(records)} contigs into {n} chunk(s)")
for path, bucket, load in zip(outputs, bins, loads):
    print(f"  {path}: {len(bucket)} contigs, {load} bp")
if records and min(loads) == 0:
    print("WARNING: at least one chunk is empty (more chunks than contigs)")
log.close()
