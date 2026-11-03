#!/usr/bin/env python3
"""
extract_domain_regions.py
VERSION: v1 | No revisions since initial creation

Extracts ONLY the HMM-matched domain region (not the whole ORF) from each
contig's best-hit protein, for either RdRp or Rep, using the coordinates
recorded by classify_domains.py's domain_classification.tsv.

Per earlier design decision: using only the conserved catalytic domain
region (rather than full ORFs, which may be much larger polyproteins with
substantial non-homologous flanking sequence) reduces alignment noise for
downstream phylogenetics, especially given how divergent these sequences
are expected to be.

Usage:
    python extract_domain_regions.py \\
        --classification domain_classification.tsv \\
        --proteins confirmed_contigs.faa \\
        --domain-type rdrp \\
        --out-fasta rdrp_domain_regions.faa
"""

import argparse
import csv
from pathlib import Path


def load_proteins(faa_path):
    """Simple FASTA parser -> dict of {protein_id: sequence}"""
    proteins = {}
    current_id = None
    current_seq = []
    with open(faa_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    proteins[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            proteins[current_id] = "".join(current_seq)
    return proteins


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", required=True)
    parser.add_argument("--proteins", required=True,
                         help="Prodigal .faa file (from domain search stage)")
    parser.add_argument("--domain-type", required=True, choices=["rdrp", "rep"])
    parser.add_argument("--out-fasta", required=True)
    args = parser.parse_args()

    proteins = load_proteins(args.proteins)

    orf_col = f"{args.domain_type}_orf_id"
    from_col = f"{args.domain_type}_ali_from"
    to_col = f"{args.domain_type}_ali_to"

    n_extracted = 0
    n_missing_protein = 0

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    with open(args.classification) as f, open(args.out_fasta, "w") as out:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            orf_id = row.get(orf_col)
            if not orf_id or orf_id == "None" or orf_id == "":
                continue  # this contig has no qualifying hit for this domain type

            ali_from = int(row[from_col])
            ali_to = int(row[to_col])

            protein_seq = proteins.get(orf_id)
            if protein_seq is None:
                n_missing_protein += 1
                print(f"[WARN] ORF '{orf_id}' listed in classification table but not "
                      f"found in {args.proteins} - skipping.")
                continue

            # ali_from/ali_to are 1-based inclusive protein coordinates from hmmsearch
            domain_seq = protein_seq[ali_from - 1:ali_to]

            contig_id = row["contig_id"]
            out.write(f">{contig_id}|{orf_id}|{args.domain_type}_domain\n")
            # wrap at 70 chars for readability, standard FASTA convention
            for i in range(0, len(domain_seq), 70):
                out.write(domain_seq[i:i+70] + "\n")
            n_extracted += 1

    print(f"Extracted {n_extracted} {args.domain_type} domain regions to {args.out_fasta}")
    if n_missing_protein:
        print(f"[WARN] {n_missing_protein} ORF ID(s) referenced in classification table "
              f"were not found in the protein FASTA - check for a mismatch between "
              f"domain search runs.")


if __name__ == "__main__":
    main()