"""
evidence_lib.py — shared parsing for the virus-ID evidence table.

Derived from Jordan's parse_evidence_table.py v4. The contig-ID normalization
here is the single most important piece: three of six tool parsers in that
pipeline silently matched ZERO contigs because each tool reports IDs in its own
format. Every parser therefore funnels its key through normalize_contig_id().

Imported by score_evidence.py and evidence_final.py via sys.path (the Snakemake
rules pass the scripts directory as params.scripts).
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

SEP = "__"   # canonical contig IDs are "{sample}{SEP}{original id}"


# ── Contig IDs ────────────────────────────────────────────────────────────────

def normalize_contig_id(raw_id):
    """Reduce any tool's contig identifier to the canonical bare ID.

    Two confirmed real-world cases:
      1. the ENTIRE original FASTA header (Cenote-Taker3 "input_name",
         VIBRANT "scaffold"):  "coasm__k141_20473 flag=0 len=2204"
      2. a call-type suffix appended by the tool (VirSorter2 "seqname"):
         "coasm__k141_61||full"
    Both rules are safe no-ops on an already-clean ID.
    """
    if not raw_id:
        return ""
    first = str(raw_id).strip().split()[0] if str(raw_id).strip() else ""
    if "||" in first:
        first = first.split("||")[0]
    return first


def sample_of(contig_id):
    """Sample (or assembly unit) a canonical contig ID belongs to."""
    return contig_id.split(SEP, 1)[0] if SEP in contig_id else ""


def load_contigs(fasta_path):
    """Ordered {contig_id: length} for one FASTA, without a Biopython import."""
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


# ── Host-domain inference from lineage strings ────────────────────────────────
# Jordan's heuristic clade tables. These map an ICTV-style lineage string to a
# likely host domain; they are approximate and deliberately conservative
# (anything unmatched stays "unknown" and is retained).

PROKARYOTE_CLADES = {
    "Caudoviricetes",      # tailed bacteriophages
    "Leviviricetes",       # ssRNA phages
    "Loebvirae",           # ssDNA phage-infecting subrealm
    "Sangervirae",         # ssDNA phage-infecting subrealm
    "Faserviricetes",      # filamentous phages
    "Vidaverviricetes",
    "Tectiliviricetes",
    "Malgrandaviricetes",  # Microviridae and relatives
}

EUKARYOTE_CLADES = {
    "Herpesvirales",
    "Nucleocytoviricota",  # NCLDV / giant eukaryotic DNA viruses
    "Pisuviricota",        # includes Picornavirales
    "Kitrinoviricota",
    "Negarnaviricota",
    "Artverviricota",
    "Duplornaviricota",
    "Cressdnaviricota",    # CRESS viruses
    "Shotokuvirae",
}


def infer_host_domain(lineage):
    """'prokaryote' | 'eukaryote' | 'unknown' for one lineage string.

    Ranks are checked most-specific first (right to left): a specific class is
    more informative than a realm.
    """
    if not lineage:
        return "unknown"
    for rank in reversed([r.strip() for r in str(lineage).split(";") if r.strip()]):
        if rank in PROKARYOTE_CLADES:
            return "prokaryote"
        if rank in EUKARYOTE_CLADES:
            return "eukaryote"
    return "unknown"


def majority_host_domain(lineages):
    """Majority vote over several lineage strings; ties and no data -> unknown."""
    votes = Counter(d for d in (infer_host_domain(x) for x in lineages) if d != "unknown")
    if not votes:
        return "unknown"
    top = votes.most_common()
    if len(top) > 1 and top[0][1] == top[1][1]:
        return "unknown"
    return top[0][0]


# ── Per-tool parsers ──────────────────────────────────────────────────────────
# Each returns {contig_id: {"called": bool, "score": float|None,
#                           "lineage": str|None}}.
# A missing or unreadable file is reported and treated as "no calls" rather
# than crashing, so one broken tool does not lose the whole run; the sanity
# checks below then fail the step if an enabled tool matched nothing.

def _rows(path, delimiter="\t"):
    if not Path(path).exists():
        print(f"[WARN] missing file: {path}", file=sys.stderr)
        return
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter=delimiter):
            yield row


def _as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_genomad(path, min_score):
    out = {}
    for row in _rows(path):
        contig = normalize_contig_id(row.get("seq_name"))
        if not contig:
            continue
        score = _as_float(row.get("virus_score") or row.get("score"))
        out[contig] = {
            "called": score is not None and score >= min_score,
            "score": score,
            "lineage": row.get("taxonomy") or None,
        }
    return out


def parse_virsorter2(path, min_score):
    """seqname carries a '||full' / '||partial' / '||lt2gene' suffix."""
    out = {}
    for row in _rows(path):
        contig = normalize_contig_id(row.get("seqname"))
        if not contig:
            continue
        score = _as_float(row.get("max_score"))
        out[contig] = {
            "called": score is not None and score >= min_score,
            "score": score,
            # max_score_group is a VirSorter2 group name, not an ICTV lineage
            "lineage": None,
            "group": row.get("max_score_group"),
        }
    return out


# VIBRANT reports a quality tier, not a probability. Mapped onto 0-1 so it can
# join the weighted score; VIBRANT lists only its positive calls.
VIBRANT_QUALITY_SCORE = {
    "complete circular": 1.0,
    "high quality draft": 0.9,
    "medium quality draft": 0.7,
    "low quality draft": 0.5,
}


def parse_vibrant(path):
    """'scaffold' holds the full original FASTA header; column is 'Quality'."""
    out = {}
    for row in _rows(path):
        raw = row.get("scaffold") or row.get("scaffold_name")
        contig = normalize_contig_id(raw)
        if not contig:
            continue
        quality = (row.get("Quality") or row.get("quality") or "").strip()
        out[contig] = {
            "called": True,          # presence == a VIBRANT viral call
            "score": VIBRANT_QUALITY_SCORE.get(quality.lower(), 0.5),
            "lineage": None,         # VIBRANT does not assign fine taxonomy
            "quality": quality,
        }
    return out


def parse_cenote_taker3(path):
    """Key on 'input_name' (full header), never 'contig' (internal ID)."""
    out = {}
    for row in _rows(path):
        raw = (row.get("input_name") or row.get("Input_name") or row.get("seq_name"))
        contig = normalize_contig_id(raw)
        if not contig:
            continue
        out[contig] = {
            "called": True,          # presence == annotated as a virus
            "score": None,           # no probability to contribute
            "lineage": row.get("taxonomy_hierarchy") or row.get("taxonomy") or None,
        }
    return out


def parse_pprmeta(path, min_score):
    out = {}
    for row in _rows(path, delimiter=","):
        contig = normalize_contig_id(row.get("Header"))
        if not contig:
            continue
        score = _as_float(row.get("Phage_score"))
        out[contig] = {
            "called": score is not None and score >= min_score,
            "score": score,
            "lineage": None,
        }
    return out


# DIAMOND / BLAST outfmt-6 column positions as requested in diamond.smk
_Q, _BITSCORE, _STAXIDS, _SCINAMES, _SKINGDOMS = 0, 11, 12, 13, 14


def parse_diamond_hits(path, top_n, bitscore_margin_pct, is_viral):
    """Majority-vote viral call per contig (Jordan's logic).

    Considers up to top_n hits, keeps those within bitscore_margin_pct of the
    contig's best bitscore, and calls viral when more than half of those
    qualifying hits are viral. `is_viral(taxid, kingdom)` decides per hit, so
    the caller supplies either a taxdump-backed test or, in test mode, the
    sskingdoms column.
    """
    out = {}
    if not Path(path).exists():
        print(f"[WARN] missing file: {path}", file=sys.stderr)
        return out

    hits = defaultdict(list)
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) <= _SKINGDOMS:
                continue
            contig = normalize_contig_id(row[_Q])
            if not contig:
                continue
            hits[contig].append({
                "bitscore": _as_float(row[_BITSCORE]) or 0.0,
                "taxids": [t for t in row[_STAXIDS].split(";") if t.strip()],
                "name": row[_SCINAMES],
                "kingdom": row[_SKINGDOMS],
            })

    for contig, contig_hits in hits.items():
        ranked = sorted(contig_hits, key=lambda h: h["bitscore"], reverse=True)[:top_n]
        if not ranked:
            continue
        cutoff = ranked[0]["bitscore"] * (bitscore_margin_pct / 100.0)
        qualifying = [h for h in ranked if h["bitscore"] >= cutoff]
        viral = [h for h in qualifying if is_viral(h["taxids"], h["kingdom"])]
        out[contig] = {
            "called": len(viral) > len(qualifying) / 2.0,
            "score": ranked[0]["bitscore"],
            "lineage": None,
            "n_qualifying_hits": len(qualifying),
            "n_viral_hits": len(viral),
            # Share of qualifying hits that look cellular: a high value on a
            # "viral" contig is the classic cellular false positive.
            "cellular_hit_fraction": round(1 - len(viral) / len(qualifying), 3) if qualifying else None,
            "top_hit": ranked[0]["name"],
            "top_taxids": ranked[0]["taxids"],
        }
    return out


# ── Sanity checks ─────────────────────────────────────────────────────────────

def check_tool_matches(tool_results, known_contigs, strict, log=None):
    """Fail when an enabled tool matched no known contig.

    This is exactly the failure Jordan hit: well-formed output, clean exit, and
    zero rows actually joined because the ID format did not match.
    """
    # Resolved at call time: a default of sys.stderr would capture the stream as
    # it was at import, before the caller redirects it to the rule's log file.
    log = log or sys.stderr
    problems = []
    for tool, results in tool_results.items():
        matched = len(set(results) & known_contigs)
        unmatched = len(set(results) - known_contigs)
        print(f"  {tool:22s} rows={len(results):6d} matched={matched:6d} unmatched={unmatched:6d}", file=log)
        if results and matched == 0:
            problems.append(f"{tool}: {len(results)} rows parsed but NONE matched a contig ID "
                            f"(ID format mismatch?)")
        elif not results:
            problems.append(f"{tool}: produced no parsable rows")
        elif unmatched > matched:
            print(f"  [WARN] {tool}: more unmatched ({unmatched}) than matched ({matched}) IDs", file=log)
    return _report(problems, strict, log)


def check_not_constant(rows, columns, strict, log=None):
    """Fail when a column is entirely one value (all False / all 0 / all empty).

    A constant column means the tool contributed nothing, which a clean exit
    code will never reveal.
    """
    log = log or sys.stderr
    problems = []
    for col in columns:
        values = {r.get(col) for r in rows}
        if len(values) <= 1 and rows:
            only = next(iter(values)) if values else None
            if only in (False, 0, 0.0, None, ""):
                problems.append(f"column '{col}' is constant ({only!r}) across all {len(rows)} contigs")
    return _report(problems, strict, log)


def _report(problems, strict, log=None):
    log = log or sys.stderr
    for p in problems:
        print(f"  [{'ERROR' if strict else 'WARN'}] {p}", file=log)
    if problems and strict:
        raise SystemExit(
            "Evidence sanity checks failed (collation.strict_checks: true):\n  - "
            + "\n  - ".join(problems)
        )
    return problems
