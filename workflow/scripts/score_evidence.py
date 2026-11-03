"""
score_evidence.py — Tier A evidence: parse every enabled detection tool into one
per-contig table, combine the evidence, and write the candidate FASTA that the
expensive Tier B steps (DIAMOND vs nr, CheckV) then run on.

Combines two approaches:
  - Jordan's union collation: a contig is a candidate if ANY tool called it,
    with confidence tiers for review priority (they never filter).
  - vine-mirus weighted scoring: a viral_score from the tools that report a
    probability, so reviewers can rank within a tier.

Called as a Snakemake script.
"""

import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.scripts)
import evidence_lib as ev   # noqa: E402

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

params = snakemake.params
tools = params.tools                  # {tool: enabled}
thresholds = params.thresholds
weights = params.weights
collation = params.collation
include_phage = bool(params.include_phage)
strict = bool(params.strict_checks)

# ── Contigs ───────────────────────────────────────────────────────────────────
lengths = {}
for path in snakemake.input.contigs:
    lengths.update(ev.load_contigs(path))
known = set(lengths)
print(f"Contig sets: {len(snakemake.input.contigs)}; contigs: {len(known)}")
if not known:
    raise SystemExit("No contigs found in the contig sets")

# ── Parse each enabled tool ───────────────────────────────────────────────────
# In test mode there is no taxdump, so viral-ness of a DIAMOND hit comes from
# the sskingdoms column. With real databases evidence_final re-derives host
# domain from the taxdump; the call here only needs viral vs not.
def _is_viral(taxids, kingdom):
    return (kingdom or "").strip().lower().startswith("virus")

dmd = params.diamond
results = {}

def _collect(tool, paths, parser):
    if not tools.get(tool):
        return
    merged = {}
    for path in paths:
        merged.update(parser(path))
    results[tool] = merged

_collect("genomad", snakemake.input.genomad,
         lambda p: ev.parse_genomad(p, thresholds["genomad_min_score"]))
_collect("virsorter2", snakemake.input.virsorter2,
         lambda p: ev.parse_virsorter2(p, thresholds["virsorter2_min_score"]))
_collect("vibrant", snakemake.input.vibrant, ev.parse_vibrant)
_collect("cenote_taker3", snakemake.input.cenote_taker3, ev.parse_cenote_taker3)
_collect("pprmeta", snakemake.input.pprmeta,
         lambda p: ev.parse_pprmeta(p, thresholds["pprmeta_min_score"]))
_collect("diamond_refseq_viral", snakemake.input.diamond_refseq_viral,
         lambda p: ev.parse_diamond_hits(p, dmd["top_n_hits"],
                                         dmd["bitscore_margin_pct"], _is_viral))

print("\nParsed tool output (rows / matched to contig IDs):")
ev.check_tool_matches(results, known, strict)

# ── Weighted score ────────────────────────────────────────────────────────────
# Only tools that report a probability carry weight. Absence from a tool's
# output is a negative call (score 0), not missing data: every enabled tool saw
# every contig, because prep_contigs applied the length gate before they ran.
scoring_tools = [t for t in ("genomad", "virsorter2", "vibrant", "pprmeta") if t in results]
weight_sum = sum(weights[t] for t in scoring_tools)
if not scoring_tools:
    print("[WARN] no score-producing tool enabled; viral_score will be 0 for every contig")
norm_weights = {t: weights[t] / weight_sum for t in scoring_tools} if weight_sum else {}
print("\nScore weights (renormalized over enabled scoring tools): "
      + ", ".join(f"{t}={w:.3f}" for t, w in norm_weights.items()))

tier_cfg = collation["confidence_tiers"]
rule = collation["candidate_rule"]
if rule not in ("union", "score", "both"):
    raise SystemExit(f"collation.candidate_rule must be union|score|both, got {rule!r}")
threshold = float(params.threshold)

rows = []
for contig, length in lengths.items():
    row = {
        "contig_id": contig,
        "sample_of_origin": ev.sample_of(contig),
        "length": length,
    }
    lineages = []
    n_called = 0
    score = 0.0
    for tool, res in results.items():
        entry = res.get(contig)
        called = bool(entry and entry.get("called"))
        row[f"{tool}_called"] = called
        row[f"{tool}_score"] = entry.get("score") if entry else None
        n_called += int(called)
        if entry and entry.get("lineage"):
            lineages.append(entry["lineage"])
        if tool in norm_weights:
            tool_score = (entry or {}).get("score")
            score += norm_weights[tool] * (tool_score if tool_score is not None else 0.0)

    row["n_tools_called"] = n_called
    row["n_tools_evaluated"] = len(results)
    row["viral_score"] = round(score, 4)

    # Tiers are review priority only.
    if n_called >= tier_cfg["tier_1_min_tools"]:
        row["confidence_tier"] = 1
    elif n_called >= tier_cfg["tier_2_min_tools"]:
        row["confidence_tier"] = 2
    elif n_called >= tier_cfg["tier_3_min_tools"]:
        row["confidence_tier"] = 3
    else:
        row["confidence_tier"] = None

    by_union = n_called >= 1
    by_score = score >= threshold
    row["is_candidate"] = {"union": by_union, "score": by_score,
                           "both": by_union and by_score}[rule]

    # Tier A host-domain signal: lineage strings from geNomad and Cenote-Taker3.
    # evidence_final refines this with the nr search.
    row["tierA_host_domain"] = ev.majority_host_domain(lineages)
    row["tierA_lineages"] = " | ".join(lineages) if lineages else None

    # Jordan's high-confidence phage rules, on Tier A evidence.
    vibrant_called = bool(row.get("vibrant_called"))
    rule_a = vibrant_called and row["tierA_host_domain"] == "prokaryote"
    rule_b = vibrant_called and n_called == 1
    row["high_confidence_phage"] = bool(rule_a or rule_b)
    row["phage_rule_triggered"] = ("rule_a_vibrant_plus_prokaryote_taxonomy" if rule_a
                                   else "rule_b_vibrant_only_call" if rule_b else None)
    rows.append(row)

# ── Sanity checks on the assembled table ──────────────────────────────────────
print("\nColumn checks:")
ev.check_not_constant(rows, [f"{t}_called" for t in results], strict)
ev.check_not_constant(rows, ["viral_score"] if norm_weights else [], strict)

# ── Write ─────────────────────────────────────────────────────────────────────
import pandas as pd   # noqa: E402  (imported late: parsing above is stdlib only)

df = pd.DataFrame(rows).sort_values(["sample_of_origin", "contig_id"])
Path(snakemake.output.scores).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(snakemake.output.scores, sep="\t", index=False)

# Candidate FASTA for Tier B. High-confidence phage is dropped here unless
# include_phage, which keeps the expensive nr/CheckV searches off obvious phage.
keep = set(df.loc[df["is_candidate"], "contig_id"])
dropped_phage = set(df.loc[df["is_candidate"] & df["high_confidence_phage"], "contig_id"])
if not include_phage:
    keep -= dropped_phage

written = 0
with open(snakemake.output.candidates, "w") as out:
    for path in snakemake.input.contigs:
        emit = False
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    emit = line[1:].strip().split()[0] in keep
                if emit:
                    out.write(line)
                    written += int(line.startswith(">"))

print(f"\nContigs evaluated        : {len(df)}")
print(f"Candidates ({rule})       : {int(df['is_candidate'].sum())}")
print(f"High-confidence phage    : {int(df['high_confidence_phage'].sum())}"
      f" ({'retained' if include_phage else 'excluded from candidates'})")
print(f"Written to candidates.fa : {written}")
print("Tier counts              : "
      + ", ".join(f"{t}={n}" for t, n in df["confidence_tier"].value_counts(dropna=False).items()))

if written == 0:
    raise SystemExit("No candidate contigs: nothing for the downstream steps to run on")
log.close()
