# Session log — vine-mirus

Append-only. Newest entry at the bottom.

## 2026-11-03

First entry. Covers the whole seven-phase merge of Jordan's "Program 1"
viral-discovery pipeline into vine-mirus, on branch `feature/program1-merge`.
Plan: `~/.claude/plans/unpack-and-read-jordan-claude-dump-10-30-ticklish-pond.md`.

### What was done

All seven planned phases, one commit each (plus two follow-ups):

| Commit | Phase |
|---|---|
| `0463c33` | 1 — globally unique contig IDs, `prep_contigs`, cleanup |
| `cd69c6c` | 2 — config restructure, preprocessing and co-assembly stages |
| `1bdc5be` | 3 — VIBRANT, Cenote-Taker3, DIAMOND searches, chunked scatter |
| `7369941` | 3 — ignore `benchmarks/` |
| `f94fe35` | 4 — evidence collation rewrite: union + weighted score + tiers |
| `f8fe713` | 5 — vOTU clustering and vOTU-level abundance |
| `fc5cb25` | 6 — database rules, legacy-CPU post-deploy, local profile |
| `f4532a7` | 7 — README, CLAUDE.md, Copilot rules, Program 2 archive |
| `d7270c4` | 7 — README DAG corrected to match the implemented rules |

Design decisions worth keeping:

- **Two tiers.** Tier A runs every enabled detection tool over all contigs,
  scattered over `n_chunks`. Tier B (DIAMOND vs nr, CheckV) runs only on
  candidates, and high-confidence phage is dropped from `candidates.fa` before
  it. nr is the single largest cost, so this gating is where the runtime budget
  goes.
- **vOTU clustering runs on `results/viral_contigs.fa`, not `candidates.fa`** —
  a deviation from the plan. On candidates, a cluster could be orphaned by its
  representative being filtered later. Consequence: the vOTU columns live in the
  summary rather than `evidence_table.tsv`, and abundance is per vOTU with
  members repeating their representative's value.
- **`votu_cluster.py` is a vendored implementation** of CheckV's
  anicalc/aniclust method. The plan assumed those scripts ship in the bioconda
  `checkv` package; checked `~/miniforge3/pkgs/checkv-1.0.1-pyhdfd78af_0/` and
  `python-scripts` is empty — they exist only in the GitHub repo. Reason is
  recorded in the script, the rule and the README.
- **Candidate rule is configurable** (`union` | `score` | `both`, default
  `union`); both numbers are always reported whichever selects.
- **Confidence tiers never filter.** They are review priority. Signals that
  could be filters are reported instead (`cellular_hit_fraction`,
  `host_domain`, `high_confidence_phage`), and unknown host domain is retained.

Runtime work delivered: DIAMOND in place of NCBI blastx, Tier B gating, phage
excluded before Tier B, minimap2 for host removal and abundance, chunked
scatter, `temp()` intermediates, vOTU representatives as the mapping reference,
and VirSorter2's nested conda env pre-built once instead of raced per chunk.
`benchmark:` directives are on the heavy rules so the gains are measurable.

### Files changed

55 files, +4450/−802 against `main`. The ones that carry the design:

- `workflow/Snakefile`: sample table, globals, `cfg_bool()`, includes, `rule all`.
- `workflow/rules/`: one `.smk` per stage; `databases.smk` holds all ten
  download/build rules plus the `databases` aggregate target.
- `workflow/scripts/evidence_lib.py`: shared parsers and
  `normalize_contig_id()` — the single defensive layer for tool ID quirks.
- `workflow/scripts/score_evidence.py` / `evidence_final.py`: Tier A then Tier B
  collation.
- `workflow/scripts/votu_cluster.py`, `run_coverm.py`, `summarize.py`.
- `workflow/scripts/program2/`: Jordan's three annotation/phylogenetics scripts
  archived unmodified (CRLF→LF only), with a README recording four open
  decisions.
- `envs/cenote_taker3.post-deploy.sh`: swaps in the SSE4.1 MMseqs2 build when
  `VINE_LEGACY_CPU=1`; as a post-deploy hook it survives env rebuilds.
- `profiles/local/config.yaml`: the 24-core box.
- `README.md`, `CLAUDE.md`, `.github/copilot-instructions.md` (last one is
  gitignored here, so it is not committed).

### Verification

- Clean-slate test-mode run from a deleted `results/`: **26/26 steps**, with
  content checks on the tables, not just the exit code.
- `snakemake --lint`: **zero** rule-level findings.
- Negative and toggle tests per phase (per-sample assembly, `pprmeta.enabled`,
  `diamond_nr.scope=all`, `votu_clustering` off, contigs-only sample sheet).
- Doc sweep: every env, script and rule file named in the docs exists.
- `git push -u origin feature/program1-merge` → pushed, tracking set.

**Not verified, and this is the important limit:** everything above is test-mode
logic. No database download rule can be exercised in test mode (tens to hundreds
of GB), so the first real run of each *is* its test. No real reads have gone
through the pipeline.

Bugs found and fixed this session, each one a case of well-formed wrong output
with a clean exit code:

- vOTU edge filter used OR across the three thresholds, which made the
  alignment-fraction threshold vacuous under `min_qcov: 0` and merged unrelated
  contigs. aniclust ANDs them. A unit test with a 250bp-shared pair exposed it:
  4 edges/2 clusters before, 2 edges/3 clusters after.
- `--config test_mode=false` silently stayed in test mode. Snakemake evaluates
  `--config` values as Python, so bare `false` arrives as the string `"false"`
  and `bool("false")` is `True` — the previously documented production command
  ran mocks. Fixed with `cfg_bool()` on every toggle. Nested keys cannot be
  `--config`-overridden at all; use `--configfile`.
- `Mean` / `Trimmed Mean` CoverM column collision: a naive `endswith()` claimed
  both columns and one sample's values overwrote another's. Fixed with
  longest-label-first matching.
- Editing a script invoked through `shell:` re-ran nothing. `evidence_lib.py`
  and `votu_cluster.py` are now declared as rule inputs with `params.scripts`
  derived from `os.path.dirname(input.lib)`.
- Config/path mismatches: `build_diamond_db` wrote `nr_viral.dmnd` while the
  config expected `refseq_viral.dmnd`; `genomad_db` pointed at the parent of
  geNomad's own nested directory.

MEGAHIT is **not** reproducible run to run: two co-assemblies of the same reads
both gave 205 contigs with different md5s, so mock draws shift and downstream
counts move by a few (144→145 viral contigs). Not a bug. The README now says to
check invariants rather than totals.

### Outstanding

- **Real-data validation on goldberg-r710** (user-run): one real sample,
  compared against Jordan's co-assembly numbers — geNomad ≈2848, VirSorter2
  ≈3357, VIBRANT ≈709, Cenote-Taker3 ≈1583. Record per-rule wall times from the
  `benchmark:` output.
- **Run each database rule once** via the `databases` target, one at a time.
  Watch the VIBRANT rule most closely: `download-db.sh` argument handling
  differs between builds, and the rule has a fallback that copies
  `databases/`/`files/` out of `$CONDA_PREFIX`.
- `profiles/slurm/config.yaml` hard-codes `slurm_partition: "main"` — check it
  against the cluster before the first submission.
- **No PR opened yet.** Branch is pushed; merging into `main` is undecided.
- `.github/copilot-instructions.md` is gitignored at `.gitignore:8`, so the
  shared source of truth does not travel with a clone. Updated locally and
  flagged rather than changing `.gitignore` unilaterally — decision pending.
- Program 2 (domain search, vConTACT3, phylogenetics) deferred by decision.
  Four open decisions in `workflow/scripts/program2/README.md`; the Rep
  reference path raises `NotImplementedError` by design.

### Rerun commands

```bash
python test/generate_test_data.py

# test-mode end-to-end (the only test there is)
snakemake --snakefile workflow/Snakefile --config samples=config/samples.tsv \
          --cores 4 --use-conda --scheduler greedy

snakemake --snakefile workflow/Snakefile --config samples=config/samples.tsv -n
snakemake --snakefile workflow/Snakefile --config samples=config/samples.tsv --lint

# databases — target name BEFORE --config
snakemake --snakefile workflow/Snakefile databases \
          --config samples=config/samples.tsv test_mode=false --cores 4 --use-conda

# production
snakemake --snakefile workflow/Snakefile --config samples=<sheet>.tsv test_mode=false \
          --profile profiles/local        # or profiles/slurm
```
