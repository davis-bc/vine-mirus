#!/usr/bin/env python3
"""
generate_test_data.py — create small synthetic paired-end FASTQ files for
pipeline testing. Produces 2 samples × ~10k read pairs each.

Sequences include a mix of:
  - Random "metagenomic" reads (background)
  - Reads derived from a synthetic phage genome (to give viral tools something to find)

Output: test/data/sample{1,2}_R{1,2}.fastq.gz, or --outdir DIR (used by
`vine-mirus test`, whose installed copy of this script cannot write beside itself).
"""

import argparse
import gzip
import random
import string
from pathlib import Path

SEED = 42
random.seed(SEED)

READ_LENGTH  = 150
N_READS      = 10_000
PHAGE_READS  = 1_000   # reads from the synthetic phage in each sample
INSERT_SIZE  = 350
INSERT_STDEV = 50
_parser = argparse.ArgumentParser(description="Write synthetic paired-end test FASTQs.")
_parser.add_argument("--outdir", type=Path, default=Path(__file__).parent / "data",
                     help="output directory (default: test/data next to this script)")
OUTDIR       = _parser.parse_args().outdir
OUTDIR.mkdir(parents=True, exist_ok=True)


# ── Helpers ───────────────────────────────────────────────────────────────────

def random_seq(length):
    return "".join(random.choices("ACGT", k=length))


def revcomp(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def make_qual(length, mean_q=35, min_q=20):
    """Return a Phred+33 quality string."""
    return "".join(
        chr(random.randint(min_q, min(40, mean_q + 5)) + 33) for _ in range(length)
    )


def reads_from_genome(genome, n, read_len, insert_size, insert_stdev):
    """Sample n paired-end reads from a circular genome."""
    pairs = []
    glen = len(genome)
    for _ in range(n):
        isize = max(read_len * 2, int(random.gauss(insert_size, insert_stdev)))
        start = random.randint(0, glen - 1)
        fragment = (genome * 3)[start: start + isize]  # wrap-around via triplication
        r1 = fragment[:read_len]
        r2 = revcomp(fragment[-read_len:])
        pairs.append((r1, r2))
    return pairs


def write_fastq_gz(path, reads, read_name_prefix):
    with gzip.open(path, "wt") as fh:
        for i, (r1, r2) in enumerate(reads):
            # For R1/R2 files we write the appropriate read from each pair
            seq = r1  # caller passes just the correct strand
            qual = make_qual(len(seq))
            fh.write(f"@{read_name_prefix}_{i+1}\n{seq}\n+\n{qual}\n")


# ── Generate synthetic phage genome (~50 kb) ─────────────────────────────────
print("Generating synthetic phage genome...")
PHAGE_GENOME_LEN = 50_000
phage_genome = random_seq(PHAGE_GENOME_LEN)

# ── Generate samples ──────────────────────────────────────────────────────────
for sample_idx in range(1, 3):
    print(f"Generating sample{sample_idx}...")

    # Background reads (random metagenome)
    bg_genome = random_seq(500_000)
    bg_pairs  = reads_from_genome(bg_genome, N_READS - PHAGE_READS, READ_LENGTH, INSERT_SIZE, INSERT_STDEV)

    # Phage-derived reads
    ph_pairs  = reads_from_genome(phage_genome, PHAGE_READS, READ_LENGTH, INSERT_SIZE, INSERT_STDEV)

    all_pairs = bg_pairs + ph_pairs
    random.shuffle(all_pairs)

    r1_reads = [(r1, r2)[0] for r1, r2 in all_pairs]
    r2_reads = [(r1, r2)[1] for r1, r2 in all_pairs]

    r1_path = OUTDIR / f"sample{sample_idx}_R1.fastq.gz"
    r2_path = OUTDIR / f"sample{sample_idx}_R2.fastq.gz"

    # Write R1
    with gzip.open(r1_path, "wt") as fh:
        for i, seq in enumerate(r1_reads):
            qual = make_qual(len(seq))
            fh.write(f"@sample{sample_idx}_read{i+1}/1\n{seq}\n+\n{qual}\n")

    # Write R2
    with gzip.open(r2_path, "wt") as fh:
        for i, seq in enumerate(r2_reads):
            qual = make_qual(len(seq))
            fh.write(f"@sample{sample_idx}_read{i+1}/2\n{seq}\n+\n{qual}\n")

    print(f"  {r1_path}  ({N_READS} read pairs)")
    print(f"  {r2_path}")

print(f"Done. Test data written to {OUTDIR}/")
