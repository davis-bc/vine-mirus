# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

@.github/copilot-instructions.md

## Running

Always pass `--snakefile workflow/Snakefile`: the Snakefile is not at the repo
root, and `configfile: "config/config.yaml"` resolves from the repo root.

```bash
# Test mode (the config default). Mocks every tool, needs no databases, seconds.
python test/generate_test_data.py   # once; test/data/*.fastq.gz is gitignored
snakemake --snakefile workflow/Snakefile --config samples=config/samples.tsv \
          --cores 4 --use-conda --scheduler greedy

# Dry run (do this before any real run)
snakemake --snakefile workflow/Snakefile --config samples=config/samples.tsv -n

# Lint — keep this at zero rule-level findings
snakemake --snakefile workflow/Snakefile --config samples=config/samples.tsv --lint

# Production, 24-core box or Slurm
snakemake --snakefile workflow/Snakefile --config samples=<sheet>.tsv test_mode=false \
          --profile profiles/local          # or --profile profiles/slurm

# Databases only (target name BEFORE --config, or Snakemake reads it as config)
snakemake --snakefile workflow/Snakefile databases \
          --config samples=config/samples.tsv test_mode=false --cores 4 --use-conda
```

There is no test suite. The test-mode run is the end-to-end check, and
"exit 0" is not the check — read the logs and the table contents. Validate one
rule by targeting its output, e.g.
`results/samples/coasm/genomad/virus_summary.tsv`.

## Architecture

Two tiers, which is where the runtime budget goes:

- **Tier A** runs every enabled detection tool over all contigs, scattered into
  `n_chunks` parallel jobs: geNomad, VirSorter2, VIBRANT, Cenote-Taker3,
  DIAMOND vs RefSeq viral, and optionally PPR-Meta. `evidence_score` collates
  them into `evidence_tierA.tsv` + `candidates.fa`.
- **Tier B** runs only on the candidates: DIAMOND vs nr and CheckV.
  `evidence_final` then assigns taxonomy, host domain and the phage rules, and
  writes `evidence_table.tsv`, `phage_evidence_table.tsv`, `viral_contigs.fa`.
- Then vOTU clustering (95% ANI / 85% AF) and CoverM against the cluster
  representatives, into `viral_contigs_abundance.tsv`.

High-confidence phage is excluded from the candidate FASTA before Tier B
(unless `phage_filtering.include_phage`), so the expensive searches skip it.

### Things that are easy to get wrong here

- **Contig IDs are namespaced**, `{unit}__{original_id}`, written by
  `prep_contigs`. MEGAHIT reuses `k141_N` in every assembly, so without this
  prefix multi-sample runs mis-join taxonomy, CheckV and abundance.
  `sample_of_origin` is derived from the prefix. `evidence_lib.normalize_contig_id`
  is the defensive layer on top: every tool reports IDs in its own format
  (VirSorter2 appends `||full`, VIBRANT and Cenote-Taker3 report the entire
  FASTA header, Cenote-Taker3's `contig` column is an unusable internal ID).
  **Every parser must funnel its key through it.** Three of six parsers in the
  pipeline this was merged from silently matched zero contigs for want of that.
- **Exit 0 is not success.** `collation.strict_checks` fails the run when an
  enabled tool matched zero contig IDs or a call column is constant; Cenote-Taker3
  gets its log marker checked because it exits 0 after doing nothing; the
  DIAMOND rule fails when every hit has `staxids=0`; the database rules verify a
  specific expected file before writing their sentinel. Preserve these when
  touching the surrounding code — each one stands for a real failure that
  produced well-formed, wrong output.
- **`cfg_bool()` in the Snakefile, not `bool()`,** for config toggles.
  Snakemake evaluates `--config` values as Python, so `test_mode=false` arrives
  as the string `"false"` and `bool("false")` is `True` — a "production" run on
  mocks. Nested keys (`votu_clustering`, `preprocessing`) cannot be overridden
  with `--config` at all; copy the config and use `--configfile`.
- **Scripts called through `shell:` are not tracked.** `evidence_lib.py` and
  `votu_cluster.py` are declared as rule inputs, with `params.scripts` derived
  from `os.path.dirname(input.lib)` so there is one source of truth (and the
  linter stays quiet). Without that, editing them re-runs nothing.
- **`temp()` on chunk files is load-bearing.** Lowering `n_chunks` once left
  stale chunks that the merge happily read, truncating 205 contigs to 50 with
  exit 0. `n_chunks` is also a rule param so a changed scatter width re-splits.
- **Java heap.** bbnorm and clumpify pre-commit their heap, so both cap `-Xmx`
  at 80% of `min(requested, MemAvailable)`. A value sized for Slurm OOM-kills a
  workstation.
- **`(zcat … || true) | head`** in `qc_reads`: under `set -o pipefail`, `head`
  closing the pipe makes `zcat` exit 141 and kills a good run.
- **`test_mode` resolves at parse time.** Any new tool rule needs all three:
  a mock command, DB inputs that become `[]`, and `envs/scripts.yaml` as its
  env in test mode (the real envs have no Python). Mocks must reproduce the
  real ID quirks, or they test nothing.

### Layout

- `workflow/rules/*.smk`, one stage each; `databases.smk` holds every download.
- `workflow/scripts/`: `evidence_lib.py` (shared parsers, imported by both
  evidence scripts), `score_evidence.py` (Tier A), `evidence_final.py` (Tier B),
  `votu_cluster.py`, `run_coverm.py`, `summarize.py`, `prep_contigs.py`,
  `split_fasta.py`, `merge_tables.py`, and the `mock_*.py` family.
- `workflow/scripts/program2/` is archived, unwired code — see its README.
- Rules reference envs as `../../envs/*.yaml` and scripts via `SCRIPTS`
  (`workflow.basedir + "/scripts"`).
- `envs/cenote_taker3.post-deploy.sh` runs on every env rebuild and swaps in the
  SSE4.1 MMseqs2 build when `VINE_LEGACY_CPU=1`. Hardware pins
  (`megahit=1.2.8`, `diamond=2.1.8` in `envs/checkv.yaml`) exist for a CPU
  without AVX2 and are harmless elsewhere — see the README's Hardware notes.

## Scientific framing

This is a screening and prioritization workflow. `evidence_table.tsv` ranks
contigs for human review and carries blank `reviewer_decision` / `reviewer_notes`
columns for that purpose. Confidence tiers are review priority and **never
filter**. Signals that could be filters are reported instead
(`cellular_hit_fraction`, `host_domain`, `high_confidence_phage`), and unknown
host domain is retained. Write about the outputs accordingly: they are not
observations that a sequence is a virus.

## Claude Code notes

- Session log: none yet. `/log-progress` will create `docs/session-log.md`.
- `.github/copilot-instructions.md` is gitignored (local-only) in this repo, so
  changes to it are not committed. Update it alongside this file anyway; flag
  the gitignore entry rather than quietly changing it.
