"""
run_taxonomy_checkv.py — extract viral contigs, assign NCBI taxonomy via DIAMOND,
and assess completeness with CheckV. Handles both real execution and test_mode mocks.
Called via Snakemake script: directive with viral_analysis conda env.
"""

import random
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

random.seed(42)

log_fh = open(snakemake.log[0], "w")
sys.stdout = log_fh
sys.stderr = log_fh

params     = snakemake.params
test_mode  = params.test_mode
threads    = snakemake.threads

# ── Extract passing viral contigs ──────────────────────────────────────────────
scores  = pd.read_csv(snakemake.input.scores, sep="\t")
passing = set(scores.loc[scores["PASS"], "contig_id"].tolist())
print(f"Passing contigs: {len(passing)}")

viral_fa = snakemake.output.viral_fa
Path(viral_fa).parent.mkdir(parents=True, exist_ok=True)

contig_to_sample = {}
with open(viral_fa, "w") as out_fa:
    for cp in snakemake.input.contigs:
        # prep_contigs writes canonical "{sample}__<id>" headers, so contig IDs
        # are globally unique and the sample is the prefix.
        for rec in SeqIO.parse(cp, "fasta"):
            if rec.id in passing and len(rec.seq) >= params.min_len:
                SeqIO.write(rec, out_fa, "fasta")
                contig_to_sample[rec.id] = rec.id.split("__", 1)[0]

print(f"Extracted {len(contig_to_sample)} viral contigs to {viral_fa}")

# ── Taxonomy ───────────────────────────────────────────────────────────────────
tax_out = snakemake.output.taxonomy
Path(tax_out).parent.mkdir(parents=True, exist_ok=True)

if test_mode:
    LINEAGES = [
        "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;Siphoviridae",
        "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;Podoviridae",
        "Viruses;Monodnaviria;Shotokuvirae;Cressdnaviricota;Repensiviricetes;;",
        "Viruses;;;;;;",
    ]
    rows = [
        {
            "contig_id": cid, "sample": s,
            "taxonomy": random.choice(LINEAGES),
            "pident": round(random.uniform(70, 100), 1),
            "evalue": float(f"1e-{random.randint(5, 50)}"),
            "bitscore": random.randint(50, 500),
        }
        for cid, s in contig_to_sample.items()
    ]
    pd.DataFrame(rows).to_csv(tax_out, sep="\t", index=False)
    print(f"[mock] Taxonomy: {len(rows)} rows")
else:
    import taxopy
    diamond_out = Path(tax_out).parent / "diamond_hits.tsv"
    subprocess.run(
        ["diamond", "blastp", "--db", params.diamond_db,
         "--query", viral_fa, "--out", str(diamond_out),
         "--outfmt", *params.outfmt.split(),
         "--evalue", str(params.evalue), "--max-target-seqs", str(params.max_target),
         "--threads", str(threads), "--more-sensitive"],
        check=True,
    )
    taxdb = taxopy.TaxDb(
        nodes_dmp=str(Path(params.taxdump) / "nodes.dmp"),
        names_dmp=str(Path(params.taxdump) / "names.dmp"),
    )
    RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    def taxid_to_lineage(taxid):
        try:
            taxon = taxopy.Taxon(int(taxid), taxdb)
            lineage = {r: "" for r in RANKS}
            for name, rank in zip(taxon.name_lineage, taxon.rank_lineage):
                if rank in lineage:
                    lineage[rank] = name
            return ";".join(lineage[r] for r in RANKS)
        except Exception:
            return ";;;;;;;"
    hits = pd.read_csv(
        diamond_out, sep="\t", header=None,
        names=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "staxids", "sscinames"],
    )
    best = hits.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
    best["taxonomy"] = best["staxids"].apply(lambda x: taxid_to_lineage(str(x).split(";")[0]))
    best["sample"]   = best["qseqid"].map(contig_to_sample)
    all_c = pd.DataFrame({"contig_id": list(contig_to_sample), "sample": list(contig_to_sample.values())})
    tax_df = all_c.merge(
        best[["qseqid", "taxonomy", "pident", "evalue", "bitscore"]].rename(columns={"qseqid": "contig_id"}),
        on="contig_id", how="left",
    )
    tax_df["taxonomy"] = tax_df["taxonomy"].fillna(";;;;;;;")
    tax_df.to_csv(tax_out, sep="\t", index=False)

# ── CheckV ─────────────────────────────────────────────────────────────────────
checkv_out = snakemake.output.checkv
checkv_dir = params.checkv_outdir
Path(checkv_dir).mkdir(parents=True, exist_ok=True)

if test_mode:
    QUALITIES = ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]
    contigs_fa = list(SeqIO.parse(viral_fa, "fasta"))
    rows = [
        {
            "contig_name": rec.id, "contig_length": len(rec.seq), "genome_copies": 1,
            "gene_count": random.randint(1, 20), "viral_genes": random.randint(1, 10),
            "host_genes": random.randint(0, 2), "checkv_quality": random.choice(QUALITIES),
            "miuvig_quality": "Genome-fragment",
            "completeness": round(random.uniform(10, 100), 1),
            "completeness_method": "AAI-based",
            "contamination": round(random.uniform(0, 5), 2),
            "provirus": "No", "proviral_length": "NA",
            "gene_count_method": "prodigal-v2.6.3",
        }
        for rec in contigs_fa
    ]
    pd.DataFrame(rows).to_csv(checkv_out, sep="\t", index=False)
    print(f"[mock] CheckV: {len(rows)} rows")
else:
    subprocess.run(
        ["checkv", "end_to_end", viral_fa, checkv_dir,
         "-d", params.checkv_db, "-t", str(threads)],
        check=True,
    )

log_fh.close()
