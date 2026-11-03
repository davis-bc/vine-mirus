#!/usr/bin/env python3
"""
classify_domains.py
VERSION: v2 | Last updated: added ali_from/ali_to coordinate tracking to
domtblout parser, for downstream domain-region extraction. v1: initial
bucket classification logic.

Parses hmmsearch --domtblout output for RdRp and Rep profile searches,
classifies each contig into one of four buckets (based on its best-hit
protein per domain type, applying the domain-completeness threshold), and
writes out both a classification table and bucket-specific FASTA subsets.

Buckets:
    rdrp_only        - RdRp domain detected (passing completeness threshold), no Rep
    rep_only         - Rep domain detected (passing completeness threshold), no RdRp
    both_detected    - BOTH detected - flagged for manual review (e.g. crucivirus-like
                        biology, or possible domain-search artifact) - routed into
                        BOTH downstream buckets
    neither_detected - Neither domain passed threshold - two-tier fallback handled
                        elsewhere (not in this script)

A contig's protein->contig mapping is inferred from Prodigal's default
header format: "<contig_id>_<orf_number> # start # end # strand # ..."

Usage:
    python classify_domains.py \\
        --rdrp-domtbl rdrp_hmmsearch.domtbl \\
        --rep-domtbl rep_hmmsearch.domtbl \\
        --contigs confirmed_contigs.fasta \\
        --min-completeness-pct 70 \\
        --out-classification domain_classification.tsv \\
        --out-rep-fasta rep_bucket_contigs.fasta \\
        --out-rdrp-fasta rdrp_bucket_contigs.fasta
"""

import argparse
import csv
from pathlib import Path
from collections import defaultdict


def prodigal_orf_to_contig(orf_id):
    """
    Prodigal default header: <contig_id>_<orf_number>
    Since contig IDs may themselves contain underscores, this assumes the
    LAST underscore-separated field is the ORF number - standard Prodigal
    convention, but worth confirming against your actual header format.
    """
    return orf_id.rsplit("_", 1)[0]


def parse_domtblout(path, min_completeness_pct):
    """
    Parses hmmsearch --domtblout format (whitespace-delimited, fixed
    columns, comment lines starting with '#').

    Relevant columns (1-indexed per HMMER docs):
        1  target name (protein/ORF ID)
        6  query name (HMM profile name)
        6  qlen (HMM profile length) - actually column 6 is tlen; qlen is column 6 in some
           versions - see note below
        16 hmm coord from (start of HMM matched)
        17 hmm coord to (end of HMM matched)

    NOTE: exact column indices for tlen/qlen can be confused between
    versions/tools - this parser uses hmmsearch's standard 23-column
    domtblout format where:
        col 1: target name, col 3: tlen (target/ORF length),
        col 4: query name, col 6: qlen (HMM profile length),
        col 16: hmm from, col 17: hmm to,
        col 18: ali from (alignment coord on the TARGET/ORF, start),
        col 19: ali to (alignment coord on the TARGET/ORF, end)
    VERIFY against your hmmer version's actual output if results look wrong.

    Returns: dict of {orf_id: {"completeness_pct": float, "ali_from": int,
    "ali_to": int}} - keeps only the single best/highest-completeness hit
    per ORF. ali_from/ali_to are the PROTEIN-coordinate positions (1-based,
    inclusive) of the matched region on the ORF's translated sequence -
    used downstream to extract just the domain region rather than the
    whole ORF.
    """
    best_hits = {}
    if not Path(path).exists():
        print(f"[WARN] domtblout file not found at {path} - treating as no hits.")
        return best_hits

    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 19:
                continue
            orf_id = fields[0]
            qlen = int(fields[5])       # HMM profile length
            hmm_from = int(fields[15])
            hmm_to = int(fields[16])
            ali_from = int(fields[17])  # protein coordinate on the ORF/target
            ali_to = int(fields[18])
            matched_len = hmm_to - hmm_from + 1
            completeness_pct = (matched_len / qlen * 100) if qlen > 0 else 0.0

            if completeness_pct >= min_completeness_pct:
                if orf_id not in best_hits or completeness_pct > best_hits[orf_id]["completeness_pct"]:
                    best_hits[orf_id] = {
                        "completeness_pct": completeness_pct,
                        "ali_from": ali_from,
                        "ali_to": ali_to,
                    }

    return best_hits


def load_all_contig_ids(fasta_path):
    ids = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def classify(rdrp_hits, rep_hits, all_contigs):
    """
    Aggregate ORF-level qualifying hits up to the contig level: a contig
    is considered to have a domain if ANY of its ORFs has a qualifying hit.
    For contigs with multiple qualifying ORFs (rare but possible on a
    multi-gene contig), the ORF with the highest completeness is selected
    as that contig's representative hit for extraction purposes.
    """
    rdrp_by_contig = defaultdict(list)
    for orf_id, hit in rdrp_hits.items():
        rdrp_by_contig[prodigal_orf_to_contig(orf_id)].append((orf_id, hit))

    rep_by_contig = defaultdict(list)
    for orf_id, hit in rep_hits.items():
        rep_by_contig[prodigal_orf_to_contig(orf_id)].append((orf_id, hit))

    def best_hit_for_contig(contig_hits):
        return max(contig_hits, key=lambda x: x[1]["completeness_pct"])

    rows = []
    for contig in all_contigs:
        has_rdrp = contig in rdrp_by_contig
        has_rep = contig in rep_by_contig

        if has_rdrp and has_rep:
            bucket = "both_detected"
        elif has_rdrp:
            bucket = "rdrp_only"
        elif has_rep:
            bucket = "rep_only"
        else:
            bucket = "neither_detected"

        row = {
            "contig_id": contig,
            "rdrp_detected": has_rdrp,
            "rep_detected": has_rep,
            "bucket": bucket,
            "rdrp_orf_id": None, "rdrp_ali_from": None, "rdrp_ali_to": None,
            "rdrp_completeness_pct": None,
            "rep_orf_id": None, "rep_ali_from": None, "rep_ali_to": None,
            "rep_completeness_pct": None,
        }

        if has_rdrp:
            orf_id, hit = best_hit_for_contig(rdrp_by_contig[contig])
            row.update({
                "rdrp_orf_id": orf_id, "rdrp_ali_from": hit["ali_from"],
                "rdrp_ali_to": hit["ali_to"], "rdrp_completeness_pct": round(hit["completeness_pct"], 1),
            })
        if has_rep:
            orf_id, hit = best_hit_for_contig(rep_by_contig[contig])
            row.update({
                "rep_orf_id": orf_id, "rep_ali_from": hit["ali_from"],
                "rep_ali_to": hit["ali_to"], "rep_completeness_pct": round(hit["completeness_pct"], 1),
            })

        rows.append(row)
    return rows


def write_fasta_subset(contigs_fasta, keep_ids, out_path):
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(contigs_fasta) as fin, open(out_path, "w") as fout:
        write_seq = False
        for line in fin:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                write_seq = seq_id in keep_ids
            if write_seq:
                fout.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rdrp-domtbl", required=True)
    parser.add_argument("--rep-domtbl", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--min-completeness-pct", type=float, default=70.0)
    parser.add_argument("--out-classification", required=True)
    parser.add_argument("--out-rep-fasta", required=True)
    parser.add_argument("--out-rdrp-fasta", required=True)
    args = parser.parse_args()

    rdrp_hits = parse_domtblout(args.rdrp_domtbl, args.min_completeness_pct)
    rep_hits = parse_domtblout(args.rep_domtbl, args.min_completeness_pct)
    all_contigs = load_all_contig_ids(args.contigs)

    rows = classify(rdrp_hits, rep_hits, all_contigs)

    Path(args.out_classification).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "contig_id", "rdrp_detected", "rep_detected", "bucket",
        "rdrp_orf_id", "rdrp_ali_from", "rdrp_ali_to", "rdrp_completeness_pct",
        "rep_orf_id", "rep_ali_from", "rep_ali_to", "rep_completeness_pct",
    ]
    with open(args.out_classification, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # rep_bucket = rep_only + both_detected (dual-domain routes into BOTH buckets)
    rep_bucket_ids = {r["contig_id"] for r in rows if r["bucket"] in ("rep_only", "both_detected")}
    rdrp_bucket_ids = {r["contig_id"] for r in rows if r["bucket"] in ("rdrp_only", "both_detected")}

    write_fasta_subset(args.contigs, rep_bucket_ids, args.out_rep_fasta)
    write_fasta_subset(args.contigs, rdrp_bucket_ids, args.out_rdrp_fasta)

    bucket_counts = defaultdict(int)
    for r in rows:
        bucket_counts[r["bucket"]] += 1

    print(f"Total contigs: {len(rows)}")
    for bucket, count in bucket_counts.items():
        print(f"  {bucket}: {count}")
    print(f"Rep bucket (rep_only + both_detected) FASTA: {len(rep_bucket_ids)} contigs "
          f"-> {args.out_rep_fasta}")
    print(f"RdRp bucket (rdrp_only + both_detected) FASTA: {len(rdrp_bucket_ids)} contigs "
          f"-> {args.out_rdrp_fasta}")


if __name__ == "__main__":
    main()