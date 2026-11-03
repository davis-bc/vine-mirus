#!/usr/bin/env python3
"""
retrieve_reference_set.py
VERSION: v2 | Last updated: implemented RdRp path (RdRp-scan). Rep path
remains deliberately stubbed (NotImplementedError) pending source
decision. v1: initial framework, both paths stubbed.

Retrieves the RdRp or Rep reference sequence set used to build supergroup-
level phylogenetic trees. This is the one step in Program 2's phylogenetics
stage that depends on internet connectivity, so it's deliberately isolated:
failure here logs to a pending-retrieval manifest and exits non-zero,
rather than being allowed to silently corrupt downstream steps or crash
the rest of the pipeline (the calling Snakefile rule uses `|| true` to
tolerate this).

SOURCE STILL TBD (see config.yaml phylogenetics.reference_retrieval):
    - RdRp: candidate sources include RdRp-scan / Serratus-associated /
      ICTV-affiliated reference alignments - not yet finalized.
    - Rep: candidate sources include published CRESS virus reference
      alignments - not yet finalized.

This script is written as a placeholder/framework: it defines the
manifest-logging and resume behavior precisely, but the actual retrieval
logic (`fetch_from_source()`) needs real implementation once a specific
source is chosen. Running this script as-is will currently always log to
the manifest and fail, by design, until that TODO is filled in.

Usage:
    python retrieve_reference_set.py \\
        --domain-type rdrp \\
        --source-config "TBD - RdRp-scan / Serratus-associated / ICTV-affiliated alignments" \\
        --out-fasta rdrp_reference_set.faa \\
        --manifest pending_reference_retrieval_manifest.tsv
"""

import argparse
import csv
import sys
from datetime import datetime
from pathlib import Path


def fetch_from_source(domain_type, source_config, out_fasta):
    """
    RDRP: implemented using RdRp-scan (Charon et al.) as the primary
    source - a curated database + HMM profiles purpose-built for divergent
    RNA virus RdRp detection/classification. Retrieved from its public
    GitHub repository. CONFIDENCE: reasonably high that this is an
    appropriate source; exact repo URL/structure should be verified at
    implementation time since it wasn't directly confirmed against live
    docs in this conversation.

    REP: NOT YET IMPLEMENTED. Unlike RdRp, there is no single well-
    established aggregator/tool serving as an obvious default source for
    CRESS/Rep reference sequences. Candidate approaches (needing a manual
    decision before this can be implemented):
        (a) Manually assemble from ICTV CRESS virus study group
            publications (Rosario/Varsani/Krupovic et al.) - curated,
            taxonomically grounded, but requires family-by-family assembly
            rather than one bulk source.
        (b) Pull RefSeq sequences matching the Rep_1 Pfam domain - broader,
            faster, but less curated/taxonomically vetted.
    This function will raise NotImplementedError for domain_type="rep"
    until one of these (or another) approach is chosen and implemented.
    """
    import subprocess

    if domain_type == "rdrp":
        Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)
        # NOTE: verify this URL/repo structure against the live RdRp-scan
        # GitHub repository before relying on this - written from general
        # knowledge of the resource, not a directly-verified current URL.
        subprocess.run(
            ["wget", "-O", out_fasta,
             "https://raw.githubusercontent.com/JustineCharon/RdRp-scan/main/RdRp_scan_db.fasta"],
            check=True
        )
        return

    raise NotImplementedError(
        f"Reference retrieval source for '{domain_type}' is not yet decided. "
        f"See function docstring for candidate approaches - this requires a "
        f"manual source decision before implementation. "
        f"Configured source note: {source_config}"
    )


def log_to_manifest(manifest_path, domain_type, source_config, error_message):
    """
    Appends a pending-retrieval entry to the manifest. If the manifest
    doesn't exist yet, creates it with a header. Designed so re-running
    this script later (once the source is implemented / connectivity is
    restored) simply succeeds and the corresponding Snakefile rule's
    output then exists, unblocking downstream alignment/tree rules -
    no upstream pipeline stages need to be rerun.
    """
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    file_exists = manifest_path.exists()

    with open(manifest_path, "a", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        if not file_exists:
            writer.writerow(["timestamp", "domain_type", "source_config", "error_message"])
        writer.writerow([
            datetime.now().isoformat(timespec="seconds"),
            domain_type,
            source_config,
            error_message,
        ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--domain-type", required=True, choices=["rdrp", "rep"])
    parser.add_argument("--source-config", required=True)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    try:
        fetch_from_source(args.domain_type, args.source_config, args.out_fasta)
        print(f"Successfully retrieved {args.domain_type} reference set -> {args.out_fasta}")
    except Exception as e:
        log_to_manifest(args.manifest, args.domain_type, args.source_config, str(e))
        print(
            f"[PENDING] Could not retrieve {args.domain_type} reference set - "
            f"logged to manifest at {args.manifest}.\n"
            f"Reason: {e}\n"
            f"This does NOT halt the rest of the pipeline - other buildable "
            f"targets will proceed. Re-run this step later once the source "
            f"is implemented/connectivity is available; downstream "
            f"alignment/tree rules for this bucket will then unblock "
            f"automatically without rerunning upstream stages.",
            file=sys.stderr
        )
        sys.exit(1)


if __name__ == "__main__":
    main()