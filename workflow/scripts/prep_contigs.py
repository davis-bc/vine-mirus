"""
prep_contigs.py — length-gate a contig FASTA and rewrite headers to canonical,
globally unique IDs. Called as a Snakemake script.

Canonical ID: "{prefix}__{first header token}", e.g. "sample1__k141_42".
The trailing MEGAHIT header text ("flag=0 multi=5.0000 len=2204") is dropped,
so every downstream tool reports the same bare ID. The prefix before "__" is
the sample (or "coasm" for co-assembly), which is how sample_of_origin is
derived downstream (see contig_origin() below / SEP in other scripts).

Fails if any ID would be duplicated after renaming, rather than silently
merging two contigs.
"""

import sys
from pathlib import Path

from Bio import SeqIO

SEP = "__"

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

prefix  = snakemake.params.prefix
min_len = int(snakemake.params.min_len)
out_fa  = snakemake.output.contigs

Path(out_fa).parent.mkdir(parents=True, exist_ok=True)

seen = set()
n_in = n_out = 0
with open(out_fa, "w") as fh:
    for rec in SeqIO.parse(snakemake.input.contigs, "fasta"):
        n_in += 1
        if len(rec.seq) < min_len:
            continue
        new_id = f"{prefix}{SEP}{rec.id}"
        if new_id in seen:
            raise ValueError(f"Duplicate contig ID after renaming: {new_id}")
        seen.add(new_id)
        fh.write(f">{new_id}\n{str(rec.seq)}\n")
        n_out += 1

print(f"{prefix}: {n_out}/{n_in} contigs >= {min_len} bp written to {out_fa}")
if n_out == 0:
    print(f"WARNING: no contigs passed the {min_len} bp gate for {prefix}")
log.close()
