#!/usr/bin/env python3
"""
mock_blastn.py — stand-in for the all-vs-all BLASTN used by vOTU clustering.

Emits BLASTN -outfmt "6 qseqid sseqid pident length qstart qend sstart send".

Besides the self-hits a real all-vs-all search produces, it makes every third
pair of consecutive contigs near-identical over almost their whole length, so
the clustering step actually collapses something and the test exercises
multi-member clusters rather than an all-singleton result.

Usage: python mock_blastn.py <contigs.fa> <out.tsv>
"""

import sys

contigs_path, out_path = sys.argv[1], sys.argv[2]

lengths = {}
cur = None
with open(contigs_path) as fh:
    for line in fh:
        if line.startswith(">"):
            cur = line[1:].strip().split()[0]
            lengths[cur] = 0
        elif cur is not None:
            lengths[cur] += len(line.strip())

ids = list(lengths)
rows = []

for contig in ids:
    length = lengths[contig]
    rows.append((contig, contig, "100.000", length, 1, length, 1, length))

# Near-identical pairs: cover 95% of each mate, which clears the default
# tcov >= 85 with room to spare.
for i in range(0, len(ids) - 1, 3):
    a, b = ids[i], ids[i + 1]
    span = int(min(lengths[a], lengths[b]) * 0.95)
    if span < 100:
        continue
    rows.append((a, b, "98.500", span, 1, span, 1, span))
    rows.append((b, a, "98.500", span, 1, span, 1, span))

with open(out_path, "w") as out:
    for row in rows:
        out.write("\t".join(str(x) for x in row) + "\n")

print(f"[mock] blastn: {len(ids)} contigs, {len(rows)} hit rows -> {out_path}",
      file=sys.stderr)
