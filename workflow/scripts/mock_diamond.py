"""
mock_diamond.py — synthetic DIAMOND blastx outfmt-6 hits for test_mode.
Called via shell: python mock_diamond.py <query_fasta> <output> <min_len>

Emits the same 15 columns the real rule requests, with REAL NCBI taxids so the
downstream majority-vote and host-domain logic is exercised: a mix of
eukaryote-infecting viruses, phage, and cellular (bacterial / fish) hits, the
last standing in for the cellular false positives nr is there to catch.
"""

import argparse
import random
from pathlib import Path

from Bio import SeqIO

random.seed(5077)

parser = argparse.ArgumentParser()
parser.add_argument("query")
parser.add_argument("output")
parser.add_argument("min_len", type=int)
args = parser.parse_args()

# (taxid, scientific name, kingdom-ish label) — real NCBI taxids
SUBJECTS = {
    "euk_virus": [
        ("12128", "Marnaviridae sp.", "Viruses"),
        ("2169561", "Genomoviridae sp.", "Viruses"),
        ("12278", "Tombusviridae sp.", "Viruses"),
    ],
    "phage": [
        ("2731619", "Caudoviricetes sp.", "Viruses"),
        ("10665", "Escherichia phage T4", "Viruses"),
    ],
    "cellular": [
        ("562", "Escherichia coli", "Bacteria"),
        ("7955", "Danio rerio", "Eukaryota"),
    ],
}

contigs = [r for r in SeqIO.parse(args.query, "fasta") if len(r.seq) >= args.min_len]
Path(args.output).parent.mkdir(parents=True, exist_ok=True)

n = 0
with open(args.output, "w") as fh:
    for rec in contigs:
        if random.random() > 0.6:          # not every contig gets a hit
            continue
        # One dominant category per contig, with occasional dissenting hits so
        # the bitscore-margin majority vote has something to resolve.
        category = random.choices(["euk_virus", "phage", "cellular"], weights=[5, 3, 2])[0]
        best = random.uniform(80, 600)
        for h in range(random.randint(1, 5)):
            cat = category if random.random() > 0.2 else random.choice(list(SUBJECTS))
            taxid, name, kingdom = random.choice(SUBJECTS[cat])
            bitscore = round(best * (1 - 0.08 * h), 1)
            fh.write("\t".join([
                rec.id,                                   # qseqid (bare contig ID)
                f"{cat}_{taxid}_{h}",                     # sseqid
                f"{random.uniform(25, 95):.1f}",          # pident
                str(random.randint(40, 400)),             # length
                str(random.randint(0, 40)),               # mismatch
                str(random.randint(0, 5)),                # gapopen
                "1", str(len(rec.seq) // 3),              # qstart qend
                "1", str(random.randint(50, 400)),        # sstart send
                f"1e-{random.randint(5, 80)}",            # evalue
                f"{bitscore}",                            # bitscore
                taxid,                                    # staxids
                name,                                     # sscinames
                kingdom,                                  # sskingdoms
            ]) + "\n")
            n += 1

print(f"[mock] DIAMOND: {n} hits for {len(contigs)} contigs -> {args.output}")
