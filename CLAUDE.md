# vine-mirus — Claude instructions

## Project overview

A Snakemake workflow for discovering viral sequences in metagenomic samples.
It takes paired-end FASTQs (or pre-assembled contigs) and produces a
reviewer-facing evidence table plus a per-vOTU abundance matrix with taxonomy.

Study context: fish tissue and lesion samples, focused on **eukaryotic**
viruses (RdRp / CRESS). Phage is excluded from the main outputs by default and
written to its own table.

## Running

Always pass `--snakefile workflow/Snakefile`: the Snakefile is not at the repo
root, and `configfile: "config/config.yaml"` resolves from the repo root.

```bash
# Test mode (the config default). Mocks every tool, needs no databases, seconds.
python test/generate_test_data.py   # once; test/data/*.fastq.gz is gitignored
snakemake --snakefile workflow/Snakefile --config samples=config/samples.tsv \
          --profile profiles/local --cores 4

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
"exit 0" is not the check — read the logs and the table contents. Passing it
means the joins and the logic work, nothing more. Validate one rule by
targeting its output, e.g. `results/samples/coasm/genomad/virus_summary.tsv`.

Per-rule resources (`threads`, `mem_mb`, `runtime`) live in the rules; the
profiles hold the executor settings and, locally, the shared `cores` / `mem_mb`
budget Snakemake packs jobs into. `profiles/slurm/config.yaml` hard-codes
`slurm_partition: "main"` — edit per cluster, and set a non-default partition
in a rule's own `resources:` when it needs one.

## Architecture

Two tiers, which is where the runtime budget goes:

```
raw FASTQs (per sample)
  └─► qc_reads: fastp → [host removal] → [clumpify dedup]   (one rule)
        └─► [bbnorm] → assemble (MEGAHIT co-assembly | per-sample | metaSPAdes)
              └─► prep_contigs: length gate + "{unit}__{id}" rename
                    ├─ TIER A, scattered over n_chunks, all contigs:
                    │    geNomad, VirSorter2, VIBRANT, Cenote-Taker3,
                    │    DIAMOND vs RefSeq viral, [PPR-Meta]
                    └─► evidence_score → evidence_tierA.tsv + candidates.fa
                          ├─ TIER B, candidates only: DIAMOND vs nr, CheckV
                          └─► evidence_final → evidence_table.tsv,
                                phage_evidence_table.tsv, viral_contigs.fa
                                └─► votu_blastn + votu_cluster (95% ANI / 85% AF)
                                      └─► coverm (minimap2-sr) on representatives
                                            └─► summary TSV
```

- **Tier A** runs every enabled detection tool over all contigs, scattered into
  `n_chunks` parallel jobs. `evidence_score` collates them into
  `evidence_tierA.tsv` + `candidates.fa`.
- **Tier B** runs only on the candidates: DIAMOND vs nr and CheckV. nr is by far
  the most expensive search, so it only sees contigs a Tier A tool already
  called. `evidence_final` then assigns taxonomy, host domain and the phage
  rules.
- Then vOTU clustering and CoverM against the cluster representatives, into
  `viral_contigs_abundance.tsv`.

High-confidence phage is excluded from the candidate FASTA before Tier B
(unless `phage_filtering.include_phage`), so the expensive searches skip it.

### Evidence collation

`score_evidence.py` (Tier A) then `evidence_final.py` (Tier B), sharing
`evidence_lib.py`.

- **Candidate rule** `collation.candidate_rule`: `union` (default — any tool
  called it), `score` (weighted `viral_score ≥ viral_score_threshold`), or
  `both`. Both numbers are always reported whichever rule selects.
- **`viral_score`** is a weighted mean over the tools that report a probability
  (geNomad `virus_score`, VirSorter2 `max_score`, VIBRANT's quality tier mapped
  to an ordinal, PPR-Meta `Phage_score`), renormalized over the enabled ones.
  Cenote-Taker3 and the DIAMOND searches contribute calls, not scores.
- **Confidence tiers** (4/2/1 tools) are **review priority only and never
  filter**.
- **Taxonomy** preference: geNomad → Cenote-Taker3 → DIAMOND best hit via the
  taxdump (taxopy); `taxonomy_source` records which was used.
- **Phage rules**: A = VIBRANT called it and the majority host domain is
  prokaryote; B = VIBRANT is the only caller.
- **Host domain** comes from lineage-clade lookup; `unknown` is retained.

### Things that are easy to get wrong here

- **Contig IDs are namespaced**, `{unit}__{original_id}` (`coasm__k141_7`,
  `sample1__k141_7`), written by `prep_contigs`, which fails on duplicates.
  MEGAHIT reuses `k141_N` in every assembly, so without this prefix
  multi-sample runs mis-join taxonomy, CheckV and abundance. Every join is on
  `contig_id`; `sample_of_origin` is derived from the prefix.
  `evidence_lib.normalize_contig_id` is the defensive layer on top: every tool
  reports IDs in its own format (VirSorter2 appends `||full`/`||partial`,
  VIBRANT's `scaffold` and Cenote-Taker3's `input_name` hold the entire FASTA
  header, Cenote-Taker3's `contig` column is an unusable internal ID).
  **Every parser must funnel its key through it.** Three of six parsers in the
  pipeline this was merged from silently matched zero contigs for want of that.
- **Exit 0 is not success.** Each check below stands for a real failure that
  produced well-formed, wrong output with a clean exit code. Preserve them when
  touching the surrounding code:
  - `collation.strict_checks` fails the run when an enabled tool matched zero
    contig IDs, or a call/score column is entirely constant.
  - Cenote-Taker3 gets its log marker checked as well as its summary file,
    because it exits 0 after doing nothing.
  - The DIAMOND rule fails when every hit has `staxids=0`.
  - `run_coverm.py` fails when CoverM's columns do not match the sample list.
  - `summarize.py` fails on contigs missing from the vOTU table or duplicated
    after a merge; `prep_contigs.py` fails on duplicate IDs.
  - Each database rule verifies a specific expected file before writing its
    `.done` sentinel.
- **Conda guard.** An `onstart:` handler in the Snakefile stops a real run
  when conda is not enabled (no `--use-conda`, no profile), because every job
  would otherwise exit 127. It skips dry-runs, `--lint` and Slurm-spawned jobs;
  `--config allow_no_conda=true` bypasses it. Keep the README's run commands
  on `--profile` or `--use-conda`. A run on atmos2 that had neither died in
  the first job with nothing more helpful than exit 127.
- **`cfg_bool()` in the Snakefile, not `bool()`,** for config toggles.
  Snakemake evaluates `--config` values as Python, so `test_mode=false` arrives
  as the string `"false"` and `bool("false")` is `True` — a "production" run on
  mocks. Nested keys (`votu_clustering`, `preprocessing`) cannot be overridden
  with `--config` at all; copy the config and use `--configfile`.
- **Scripts called through `shell:` are not tracked.** `evidence_lib.py` and
  `votu_cluster.py` are declared as rule inputs, with `params.scripts` derived
  from `os.path.dirname(input.lib)` so there is one source of truth (and the
  linter stays quiet). Without that, editing them re-runs nothing.
- **Derive path-prefix params from inputs/outputs**
  (`lambda wc, output: os.path.dirname(output.done)`), never hardcode them —
  the linter flags it, and it breaks without a shared filesystem.
- **`temp()` on chunk files is load-bearing.** Lowering `n_chunks` once left
  stale chunks that the merge happily read, truncating 205 contigs to 50 with
  exit 0. `n_chunks` is also a rule param so a changed scatter width re-splits.
- **Java heap.** bbnorm and clumpify pre-commit their heap, so both cap `-Xmx`
  at 80% of `min(requested, MemAvailable)`. A value sized for Slurm OOM-kills a
  workstation.
- **`(zcat … || true) | head`** in `qc_reads`: under `set -o pipefail`, `head`
  closing the pipe makes `zcat` exit 141 and kills a good run.
- **`test_mode` resolves at parse time.** Any new tool rule needs all three:
  a mock command selected by the module-level `TEST_MODE` conditional, DB
  sentinel inputs that become `[]`, and `envs/scripts.yaml` as its env in test
  mode (the real envs have no Python). Mocks must reproduce the real ID quirks,
  or they test nothing.

### Databases

`databases.smk` owns all of them, under `databases/`, each gated by a `.done`
sentinel that tool rules declare as a conditional input.

Several downloaders write one level deeper than the consuming tool looks:
`genomad download-database X` → `X/genomad_db/`; `checkv download_database X` →
`X/checkv-db-v1.N/` (resolved at run time); `download-db.sh X` → `X/databases/`
plus a sibling `files/` from the package; `virsorter setup` ignores `--db-dir`
and writes to the **current directory**, so its rule `cd`s into the target.

DIAMOND databases are built with `--taxonmap/--taxonnodes/--taxonnames` so hits
carry real taxids. nr (~250 GB) is never fetched implicitly: without
`download_nr: true` a rule fails explaining the three ways to provide it.

### Conda environments

One env per genuine dependency island (more than the ideal 4–5 because the
conflicts are real: CheckV caps DIAMOND at ≤2.1.8 while the searches want a
current release; VIBRANT needs `numpy<1.24`; PPR-Meta needs TensorFlow). Check
whether a tool fits an existing env before adding one. In test mode every tool
rule uses `envs/scripts.yaml`.

Hardware: the reference box (`goldberg-r710`, 24 cores) has SSE4.2 but no
AVX2, where several current conda builds die with "Illegal instruction".
`megahit=1.2.8` and `diamond=2.1.8` (in `envs/checkv.yaml`) are pinned
unconditionally, harmless on newer CPUs. MMseqs2 has no suitable conda build,
so `envs/cenote_taker3.post-deploy.sh` swaps in the static SSE4.1 binary when
`VINE_LEGACY_CPU=1`; as a post-deploy hook it re-applies on every env rebuild,
which a manual swap does not.

### Layout

- `workflow/Snakefile`: sample table, globals, `cfg_bool()`, includes, `rule all`.
- `workflow/rules/*.smk`, one stage each; `databases.smk` holds every download
  and the `databases` target.
- `workflow/scripts/`: `evidence_lib.py` (shared parsers, imported by both
  evidence scripts), `score_evidence.py` (Tier A), `evidence_final.py` (Tier B),
  `votu_cluster.py`, `run_coverm.py`, `summarize.py`, `prep_contigs.py`,
  `split_fasta.py`, `merge_tables.py`, and the `mock_*.py` family.
- `workflow/scripts/program2/` is archived, unwired code — see its README.
- Rules reference envs as `../../envs/*.yaml` and scripts via `SCRIPTS`
  (`workflow.basedir + "/scripts"`).
- `config/` (`config.yaml`, `samples.tsv`), `profiles/{slurm,local}/`,
  `test/` (`generate_test_data.py`, gitignored `data/`). `results/`, `logs/`,
  `databases/`, `benchmarks/` are all gitignored.

### Outputs

Terminal: `results/evidence_table.tsv` (the reviewer table — every contig
screened), `results/phage_evidence_table.tsv`, `results/viral_contigs.fa`,
`results/viral_contigs_abundance.tsv`, `results/qc/multiqc_report.html`.

The summary carries `contig_id`, `sample_of_origin`, `length`, `viral_score`,
`n_tools_called`, `confidence_tier`, `taxonomy`, `host_domain`,
`checkv_quality`, `checkv_completeness`, the vOTU columns, then one abundance
column per sample. Abundance is **per vOTU**: members repeat their
representative's value, and `is_votu_representative` marks which contig carried
the mapping.

## Scientific framing

This is a screening and prioritization workflow. `evidence_table.tsv` ranks
contigs for human review and carries blank `reviewer_decision` / `reviewer_notes`
columns for that purpose. Confidence tiers are review priority and **never
filter**. Signals that could be filters are reported instead
(`cellular_hit_fraction`, `host_domain`, `high_confidence_phage`), and unknown
host domain is retained. Write code and prose about the outputs accordingly:
they are not observations that a sequence is a virus.

## Claude Code notes

- Session log: `docs/session-log.md` (single-file convention, appended by
  `/log-progress`).
