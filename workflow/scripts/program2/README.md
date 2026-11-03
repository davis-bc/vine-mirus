# Program 2 scripts (archived, not wired into the workflow)

These three scripts come from Jordan's "Program 2" — the annotation and
phylogenetics stage that follows Program 1's discovery screen. They are kept
here **unmodified** (apart from converting CRLF line endings to LF) so the work
is not lost, and because the open design decisions below are easier to settle
with the code in front of you.

No Snakemake rule references them. Nothing in `workflow/Snakefile` will run
them, and they are not covered by the test-mode run.

| Script | What it does |
|--------|--------------|
| `classify_domains.py` | Parses `hmmsearch --domtblout` for RdRp and Rep profiles and sorts each contig into one of four buckets: `rdrp_only`, `rep_only`, `both_detected` (flagged for manual review — crucivirus-like biology, or a search artifact), `neither_detected`. Also writes bucket FASTA subsets and records `ali_from`/`ali_to`. |
| `extract_domain_regions.py` | Pulls out only the HMM-matched domain region of each best-hit protein, using the coordinates `classify_domains.py` recorded. Aligning the conserved catalytic domain rather than whole ORFs cuts noise, which matters for sequences this divergent. |
| `retrieve_reference_set.py` | Fetches the RdRp or Rep reference set used to build supergroup-level trees. The RdRp path is implemented against RdRp-scan; **the Rep path raises `NotImplementedError` by design**, pending the source decision below. It is the one internet-dependent step, deliberately isolated so a failure logs to a manifest instead of corrupting downstream steps. |

## Open decisions, to settle before wiring any of this up

1. **Rep reference source.** Undecided. Candidates are published CRESS virus
   reference alignments. `retrieve_reference_set.py` will always fail on the Rep
   path until a source is chosen and `fetch_from_source()` is implemented. The
   RdRp source (RdRp-scan / Serratus-associated / ICTV-affiliated alignments)
   is also marked "not yet finalized" in the script's own header.
2. **Domain-completeness threshold (70%).** A placeholder. It decides which
   contigs reach the trees at all, so it wants benchmarking on real data before
   it is treated as a setting rather than a guess.
3. **Clade-size threshold.** The minimum clade size for a supergroup tree was
   never fixed.
4. **vConTACT3 scope.** Jordan's config runs it on the `rep_only` and
   `both_detected` buckets only. The reasoning, from the vConTACT3 paper
   (Bolduc et al. 2025): genus-level accuracy is strong for
   eukaryote-infecting Monodnaviria/CRESS (96.3%), but genus through order are
   all "not possible" for eukaryote-infecting Riboviria — the realm holding most
   of this project's RNA virus content. RdRp-bucket contigs therefore rely on
   phylogenetics instead, which is why that stage exists.

## How these would attach

Program 1 ends at `results/evidence_table.tsv` plus `results/viral_contigs.fa`
and the per-vOTU abundance table. Program 2 starts from the contigs a human has
confirmed in the reviewer table — so the natural entry point is a new sample
sheet or a filter on `reviewer_decision`, not an automatic hand-off. Keeping the
manual review in the middle is deliberate: the screen produces candidates for
assessment, not conclusions.
