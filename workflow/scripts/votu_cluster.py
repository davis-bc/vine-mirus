#!/usr/bin/env python3
"""
votu_cluster.py — dereplicate viral contigs into vOTUs (viral operational
taxonomic units) by ANI / alignment-fraction greedy centroid clustering.

This implements the method CheckV documents in its `anicalc.py` + `aniclust.py`
utilities, which is the community convention for vOTU definition (95% ANI over
85% of the shorter sequence; Gregory et al., Global Ocean Virome). Those two
scripts live in the CheckV *repository*, not in the bioconda `checkv` package,
so they cannot be invoked from the conda env the pipeline builds — hence this
vendored implementation. It follows their semantics deliberately:

  per ordered pair (query, target), from all-vs-all BLASTN HSPs:
    ani  = length-weighted mean percent identity over all HSPs
    qcov = % of the QUERY  covered by the merged HSP intervals
    tcov = % of the TARGET covered by the merged HSP intervals
  an edge query -> target is kept when ALL of
    ani >= min_ani, tcov >= min_tcov, qcov >= min_qcov
  The conditions are ANDed, as aniclust does; min_qcov defaults to 0, which
  means "impose no query-coverage requirement". Treating them as alternatives
  instead would make min_qcov: 0 satisfy the filter on its own and void the
  alignment-fraction threshold entirely.
  clustering is greedy and longest-first: each still-unassigned contig becomes
  a representative and absorbs every unassigned target it has an edge to.

Intervals are merged per sequence before coverage is computed, so overlapping
HSPs cannot inflate the aligned length past the sequence itself. Identity is
weighted by each HSP's own (unmerged) length, as anicalc does; this is an
approximation when HSPs overlap heavily.

Usage (called from workflow/rules/votu.smk):
    python votu_cluster.py --blast all_vs_all.tsv --contigs viral_contigs.fa \
        --min-ani 95 --min-tcov 85 --min-qcov 0 \
        --out-clusters votu_clusters.tsv --out-reps votu_representatives.fa
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

# BLASTN -outfmt "6 qseqid sseqid pident length qstart qend sstart send"
N_FIELDS = 8


def load_contigs(fasta_path):
    """Ordered {contig_id: length}."""
    lengths = {}
    cur = None
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].strip().split()[0]
                lengths[cur] = 0
            elif cur is not None:
                lengths[cur] += len(line.strip())
    return lengths


def merge_intervals(intervals):
    """Merge overlapping/abutting 1-based inclusive (start, end) intervals."""
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def _span(merged):
    return sum(end - start + 1 for start, end in merged)


def parse_blast(path, known):
    """{(query, target): {'q': [...], 't': [...], 'pid_sum': x, 'len_sum': y}}.

    Self-hits are skipped. Hits naming a contig absent from the FASTA are
    counted and reported rather than silently dropped.
    """
    pairs = defaultdict(lambda: {"q": [], "t": [], "pid_sum": 0.0, "len_sum": 0})
    n_rows = n_self = n_unknown = 0
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < N_FIELDS:
                continue
            n_rows += 1
            query, target = row[0], row[1]
            if query == target:
                n_self += 1
                continue
            if query not in known or target not in known:
                n_unknown += 1
                continue
            try:
                pident, length = float(row[2]), int(row[3])
                qstart, qend, tstart, tend = (int(x) for x in row[4:8])
            except ValueError:
                continue
            entry = pairs[(query, target)]
            entry["q"].append((min(qstart, qend), max(qstart, qend)))
            entry["t"].append((min(tstart, tend), max(tstart, tend)))
            entry["pid_sum"] += pident * length
            entry["len_sum"] += length
    return pairs, n_rows, n_self, n_unknown


def pair_metrics(pairs, lengths):
    """{(query, target): (ani, qcov, tcov)} with coverages as percentages."""
    out = {}
    for (query, target), entry in pairs.items():
        if not entry["len_sum"]:
            continue
        ani = entry["pid_sum"] / entry["len_sum"]
        qlen, tlen = lengths[query], lengths[target]
        qcov = 100.0 * _span(merge_intervals(entry["q"])) / qlen if qlen else 0.0
        tcov = 100.0 * _span(merge_intervals(entry["t"])) / tlen if tlen else 0.0
        # Merged coverage cannot exceed the sequence; clamp float drift only.
        out[(query, target)] = (ani, min(qcov, 100.0), min(tcov, 100.0))
    return out


def build_edges(metrics, min_ani, min_tcov, min_qcov):
    edges = defaultdict(list)
    for (query, target), (ani, qcov, tcov) in metrics.items():
        if ani >= min_ani and tcov >= min_tcov and qcov >= min_qcov:
            edges[query].append(target)
    return edges


def greedy_cluster(lengths, edges):
    """Longest-first greedy centroid clustering -> {representative: [members]}.

    Ties on length break on contig ID so the result is reproducible.
    """
    order = sorted(lengths, key=lambda c: (-lengths[c], c))
    assigned = set()
    clusters = {}
    for contig in order:
        if contig in assigned:
            continue
        members = [contig]
        assigned.add(contig)
        for target in edges.get(contig, ()):
            if target not in assigned:
                members.append(target)
                assigned.add(target)
        clusters[contig] = members
    return clusters


def write_clusters(clusters, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig_id", "votu_representative", "votu_cluster_size",
                         "is_votu_representative"])
        for rep, members in clusters.items():
            for member in members:
                writer.writerow([member, rep, len(members), member == rep])


def write_reps(clusters, contigs_fasta, path):
    keep = set(clusters)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    written = 0
    emit = False
    with open(contigs_fasta) as fin, open(path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                emit = line[1:].strip().split()[0] in keep
                written += int(emit)
            if emit:
                fout.write(line)
    return written


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--min-ani", type=float, default=95.0)
    parser.add_argument("--min-tcov", type=float, default=85.0)
    parser.add_argument("--min-qcov", type=float, default=0.0)
    parser.add_argument("--out-clusters", required=True)
    parser.add_argument("--out-reps", required=True)
    args = parser.parse_args()

    lengths = load_contigs(args.contigs)
    if not lengths:
        sys.exit(f"No contigs in {args.contigs}: nothing to cluster")

    pairs, n_rows, n_self, n_unknown = parse_blast(args.blast, set(lengths))
    print(f"Contigs            : {len(lengths)}")
    print(f"BLASTN rows        : {n_rows} ({n_self} self-hits skipped)")
    if n_unknown:
        print(f"[WARN] {n_unknown} BLASTN rows named a contig absent from {args.contigs}")
    print(f"Contig pairs        : {len(pairs)}")
    if n_rows and not pairs:
        # Every contig is its own vOTU. Legitimate on a good co-assembly, but
        # also what an ID-format mismatch between the FASTA and the BLAST
        # output looks like, so say so explicitly.
        print("[WARN] no between-contig alignments: every contig will be a singleton vOTU")

    metrics = pair_metrics(pairs, lengths)
    edges = build_edges(metrics, args.min_ani, args.min_tcov, args.min_qcov)
    print(f"Edges kept (ani>={args.min_ani}, tcov>={args.min_tcov}, qcov>={args.min_qcov}): "
          f"{sum(len(v) for v in edges.values())}")

    clusters = greedy_cluster(lengths, edges)
    n_members = sum(len(m) for m in clusters.values())
    if n_members != len(lengths):
        sys.exit(f"Clustering lost contigs: {n_members} members for {len(lengths)} contigs")

    write_clusters(clusters, args.out_clusters)
    written = write_reps(clusters, args.contigs, args.out_reps)
    if written != len(clusters):
        sys.exit(f"Wrote {written} representative sequences for {len(clusters)} clusters")

    singletons = sum(1 for m in clusters.values() if len(m) == 1)
    largest = max((len(m) for m in clusters.values()), default=0)
    print(f"\nvOTUs (clusters)   : {len(clusters)}")
    print(f"Singletons         : {singletons} ({100.0 * singletons / len(clusters):.1f}%)")
    print(f"Largest cluster    : {largest} contigs")
    print(f"Dereplication      : {len(lengths)} -> {len(clusters)} "
          f"({100.0 * (1 - len(clusters) / len(lengths)):.1f}% collapsed)")
    print(f"Clusters -> {args.out_clusters}")
    print(f"Representatives -> {args.out_reps}")


if __name__ == "__main__":
    main()
