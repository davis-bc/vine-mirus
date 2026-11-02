"""
evidence_final.py — Tier B: add the nr search and CheckV to the Tier A table,
assign a taxonomic lineage, finalize the host-domain and phage calls, and write
the reviewer-facing tables plus the viral contig FASTA.

Taxonomy preference order (most to least specific source):
  1. geNomad lineage   - ICTV ranks, assigned by a viral classifier
  2. Cenote-Taker3     - taxonomy_hierarchy
  3. DIAMOND best hit  - lineage resolved from the hit's taxid via the taxdump
Writes which source was used (taxonomy_source) so the provenance is visible.

Called as a Snakemake script.
"""

import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.scripts)
import evidence_lib as ev   # noqa: E402

import pandas as pd   # noqa: E402

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

params = snakemake.params
test_mode = bool(params.test_mode)
strict = bool(params.strict_checks)
include_phage = bool(params.include_phage)
dmd = params.diamond

df = pd.read_csv(snakemake.input.scores, sep="\t")
print(f"Tier A table: {len(df)} contigs, {int(df['is_candidate'].sum())} candidates")

# ── Taxonomy backend ──────────────────────────────────────────────────────────
# With real databases, lineages come from the NCBI taxdump via taxopy. In test
# mode there is no taxdump, so the DIAMOND columns (sscinames / sskingdoms)
# stand in and only the viral/cellular distinction is available.
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
EMPTY_LINEAGE = ";" * (len(RANKS) - 1)

if test_mode or not params.taxdump:
    taxdb = None
    print("[test mode] no taxdump; DIAMOND taxonomy limited to the sskingdoms column")
else:
    import taxopy
    taxdb = taxopy.TaxDb(nodes_dmp=str(Path(params.taxdump) / "nodes.dmp"),
                         names_dmp=str(Path(params.taxdump) / "names.dmp"))


def taxid_lineage(taxid):
    """';'-joined lineage at the standard ranks, or empty when unresolvable."""
    if taxdb is None or not taxid:
        return ""
    try:
        taxon = taxopy.Taxon(int(taxid), taxdb)
    except Exception:
        return ""
    ranks = {r: "" for r in RANKS}
    for name, rank in zip(taxon.name_lineage, taxon.rank_lineage):
        if rank in ranks:
            ranks[rank] = name
    return ";".join(ranks[r] for r in RANKS)


def hit_is_viral(taxids, kingdom):
    if taxdb is not None:
        for taxid in taxids:
            try:
                if "10239" in taxopy.Taxon(int(taxid), taxdb).taxid_lineage_str.split(";"):
                    return True
            except Exception:
                continue
        return False
    return (kingdom or "").strip().lower().startswith("virus")


# ── nr search (candidates only by default) ────────────────────────────────────
nr = {}
if snakemake.input.diamond_nr:
    for path in snakemake.input.diamond_nr:
        nr.update(ev.parse_diamond_hits(path, dmd["top_n_hits"],
                                        dmd["bitscore_margin_pct"], hit_is_viral))
    candidates = set(df.loc[df["is_candidate"], "contig_id"])
    print(f"\nnr hits parsed for {len(nr)} contigs "
          f"({len(set(nr) & candidates)} of {len(candidates)} candidates)")
    ev.check_tool_matches({"diamond_nr": nr}, set(df["contig_id"]), strict)
else:
    print("\nnr search disabled")

df["diamond_nr_called"] = df["contig_id"].map(lambda c: bool(nr.get(c, {}).get("called")))
df["diamond_nr_top_hit"] = df["contig_id"].map(lambda c: nr.get(c, {}).get("top_hit"))
# High cellular_hit_fraction on a called contig is the classic cellular false
# positive: reported, never auto-filtered.
df["cellular_hit_fraction"] = df["contig_id"].map(
    lambda c: nr.get(c, {}).get("cellular_hit_fraction"))

if nr:
    # nr is an independent line of evidence, so it joins the tool count/tier.
    df["n_tools_called"] = df["n_tools_called"] + df["diamond_nr_called"].astype(int)
    df["n_tools_evaluated"] = df["n_tools_evaluated"] + 1
    tiers = params.confidence_tiers

    def tier_of(n):
        if n >= tiers["tier_1_min_tools"]:
            return 1
        if n >= tiers["tier_2_min_tools"]:
            return 2
        if n >= tiers["tier_3_min_tools"]:
            return 3
        return None

    df["confidence_tier"] = df["n_tools_called"].map(tier_of)

# ── CheckV ────────────────────────────────────────────────────────────────────
checkv = pd.read_csv(snakemake.input.checkv, sep="\t")
checkv["contig_id"] = checkv["contig_id"].map(ev.normalize_contig_id)
keep_cols = [c for c in ["contig_id", "checkv_quality", "completeness",
                         "contamination", "provirus"] if c in checkv.columns]
checkv = checkv[keep_cols].rename(columns={"completeness": "checkv_completeness",
                                           "contamination": "checkv_contamination",
                                           "provirus": "checkv_provirus"})
matched = len(set(checkv["contig_id"]) & set(df["contig_id"]))
print(f"\nCheckV: {len(checkv)} rows, {matched} matched to contig IDs")
if len(checkv) and matched == 0:
    ev._report(["CheckV rows matched no contig ID"], strict)
df = df.merge(checkv, on="contig_id", how="left")

# ── Final taxonomy and host domain ────────────────────────────────────────────
def pick_taxonomy(row):
    """geNomad > Cenote-Taker3 > DIAMOND, using Tier A lineages where present."""
    raw = row.get("tierA_lineages")
    if isinstance(raw, str) and raw.strip():
        first = raw.split(" | ")[0].strip()
        if first:
            source = "genomad" if row.get("genomad_called") else "cenote_taker3"
            return first, source
    nr_entry = nr.get(row["contig_id"], {})
    if nr_entry:
        lineage = taxid_lineage((nr_entry.get("top_taxids") or [None])[0]) if taxdb else ""
        if lineage:
            return lineage, "diamond_nr"
        if nr_entry.get("top_hit"):
            return str(nr_entry["top_hit"]), "diamond_nr_top_hit_name"
    return EMPTY_LINEAGE, "unassigned"


tax = df.apply(pick_taxonomy, axis=1, result_type="expand")
df["taxonomy"], df["taxonomy_source"] = tax[0], tax[1]

def final_domain(row):
    """Majority vote over every lineage available for this contig."""
    lineages = []
    if isinstance(row.get("tierA_lineages"), str):
        lineages += [x for x in row["tierA_lineages"].split(" | ") if x]
    if isinstance(row.get("taxonomy"), str) and row.get("taxonomy_source") == "diamond_nr":
        lineages.append(row["taxonomy"])
    return ev.majority_host_domain(lineages)


df["host_domain"] = df.apply(final_domain, axis=1)

# Jordan's rules, re-evaluated with the refined host domain.
vibrant_called = df["vibrant_called"] if "vibrant_called" in df.columns else False
rule_a = vibrant_called & (df["host_domain"] == "prokaryote")
rule_b = vibrant_called & (df["n_tools_called"] == 1)
df["high_confidence_phage"] = rule_a | rule_b
df["phage_rule_triggered"] = None
df.loc[rule_b, "phage_rule_triggered"] = "rule_b_vibrant_only_call"
df.loc[rule_a, "phage_rule_triggered"] = "rule_a_vibrant_plus_prokaryote_taxonomy"

df["host_domain_classification"] = "unknown_retained_by_default"
df.loc[df["host_domain"] == "eukaryote", "host_domain_classification"] = "eukaryote_likely"
df.loc[df["high_confidence_phage"], "host_domain_classification"] = "prokaryote_likely_phage"

# Reviewer columns, filled in by hand after the run.
df["reviewer_decision"] = ""
df["reviewer_notes"] = ""

# ── Write ─────────────────────────────────────────────────────────────────────
lead = ["contig_id", "sample_of_origin", "length", "viral_score", "n_tools_called",
        "n_tools_evaluated", "confidence_tier", "is_candidate", "taxonomy",
        "taxonomy_source", "host_domain", "host_domain_classification",
        "checkv_quality", "checkv_completeness", "checkv_contamination",
        "high_confidence_phage", "phage_rule_triggered", "cellular_hit_fraction"]
ordered = [c for c in lead if c in df.columns] + [c for c in df.columns if c not in lead]
df = df[ordered].sort_values(["confidence_tier", "viral_score"], ascending=[True, False])

for path in (snakemake.output.evidence, snakemake.output.phage):
    Path(path).parent.mkdir(parents=True, exist_ok=True)

df.to_csv(snakemake.output.evidence, sep="\t", index=False)
df[df["high_confidence_phage"]].to_csv(snakemake.output.phage, sep="\t", index=False)

keep = set(df.loc[df["is_candidate"], "contig_id"])
if not include_phage:
    keep -= set(df.loc[df["high_confidence_phage"], "contig_id"])

written = 0
emit = False
with open(snakemake.output.viral_fa, "w") as out, open(snakemake.input.candidates) as fh:
    for line in fh:
        if line.startswith(">"):
            emit = line[1:].strip().split()[0] in keep
        if emit:
            out.write(line)
            written += int(line.startswith(">"))

print(f"\nEvidence table      : {len(df)} rows -> {snakemake.output.evidence}")
print(f"High-confidence phage: {int(df['high_confidence_phage'].sum())} -> {snakemake.output.phage}")
print(f"Viral contigs        : {written} -> {snakemake.output.viral_fa}")
print("Taxonomy sources     : "
      + ", ".join(f"{k}={v}" for k, v in df["taxonomy_source"].value_counts().items()))
print("Host domain          : "
      + ", ".join(f"{k}={v}" for k, v in df["host_domain"].value_counts().items()))

if written == 0:
    raise SystemExit("No viral contigs written: downstream abundance would be empty")
log.close()
