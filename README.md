# vine-mirus

**Viral metagenome discovery workflow** — an end-to-end Snakemake pipeline for
identifying, annotating, and quantifying viral contigs in metagenomic samples.

---

## Overview

vine-mirus assembles paired-end reads (or accepts pre-assembled contigs), screens
them through three complementary viral detection tools, fuses their scores into a
single evidence-weighted probability, assigns NCBI taxonomy, assesses genome
completeness, and reports per-contig abundances across all samples.

```
FASTQs (or contigs)
       │
       ▼
fastp (QC)  ──►  MEGAHIT (assembly)
                        │
          ┌─────────────┼─────────────┐
          ▼             ▼             ▼
   geNomad  VirSorter2  VIBRANT  Cenote-Taker3  DIAMOND/RefSeq  [PPR-Meta]
          └─────────────┴──────┬──────┴───────────────┘
                        ▼
               evidence_score  ←  Tier A: union + weighted score + tiers
                        │                 → candidates.fa
                        ▼
          DIAMOND vs nr  +  CheckV  ←  Tier B: candidates only
                        │
                        ▼
               evidence_final  ←  taxonomy, host domain, phage rules
                        │
                        ▼
                      vOTU  ←  95% ANI / 85% AF → representatives
                        │
                        ▼
                     CoverM  ←  minimap2-sr → TPM per vOTU per sample
                        │
                        ▼
          results/viral_contigs_abundance.tsv
          results/evidence_table.tsv
          results/viral_contigs.fa
```

**Databases are downloaded automatically** on first execution and cached in
`databases/` — no manual setup required.

---

## Requirements

| Tool | Minimum version |
|------|----------------|
| [Conda / Mamba](https://github.com/conda-forge/miniforge) | any recent |
| [Snakemake](https://snakemake.readthedocs.io) | ≥ 9.0 |
| snakemake-executor-plugin-slurm | ≥ 0.4 (SLURM only) |

Install a Snakemake environment with conda:

```bash
mamba create -n snakemake -c conda-forge -c bioconda \
    snakemake snakemake-executor-plugin-slurm
conda activate snakemake
```

---

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/davis-bc/vine-mirus.git
cd vine-mirus

# 2. Create your sample sheet
cp config/samples.tsv config/my_samples.tsv
# Edit my_samples.tsv — see "Sample Sheet" section below

# 3. Run (databases download automatically on first run)
snakemake --snakefile workflow/Snakefile \
          --config samples=config/my_samples.tsv \
          --cores 16 --use-conda

# 4. Check results
ls results/
# results/viral_contigs_abundance.tsv  ← main output
# results/viral_contigs.fa             ← all passing viral sequences
```

---

## Sample Sheet

The sample sheet is a TSV with a `sample` column and one of the following input
format options.

### Option 1 — Paired-end reads (full pipeline)

```tsv
sample  r1                              r2
sample1 /data/sample1_R1.fastq.gz       /data/sample1_R2.fastq.gz
sample2 /data/sample2_R1.fastq.gz       /data/sample2_R2.fastq.gz
```

Runs: QC (fastp → optional host removal → optional dedup) → normalization → co-assembly (or per-sample assembly) → viral detection → taxonomy/CheckV → CoverM.

In `coassembly` mode, contig IDs are prefixed `coasm__`. In `per_sample` mode they are prefixed `{sample}__`.

### Option 2 — Pre-assembled contigs (viral detection only)

```tsv
sample  contigs
sample1 /data/sample1_contigs.fa
sample2 /data/sample2_contigs.fa
```

Skips QC and assembly. Viral detection, taxonomy, and CheckV are run but no
read-based abundance estimation is performed (CoverM is skipped).

### Option 3 — Pre-assembled contigs + reads for abundance

```tsv
sample  r1                              r2                              contigs
sample1 /data/sample1_R1.fastq.gz       /data/sample1_R2.fastq.gz       /data/sample1.fa
sample2 /data/sample2_R1.fastq.gz       /data/sample2_R2.fastq.gz       /data/sample2.fa
```

Uses provided contigs for viral detection but still runs fastp + CoverM for
read-level abundance across all samples.

---

## Configuration

All parameters are in `config/config.yaml`. Key settings:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_contig_length` | `1000` | Minimum contig length passed to viral tools |
| `viral_score_threshold` | `0.5` | Minimum weighted evidence score to retain a contig |
| `tool_weights.genomad` | `0.34` | Weight for geNomad score |
| `tool_weights.virsorter2` | `0.33` | Weight for VirSorter2 score |
| `tool_weights.pprmeta` | `0.33` | Weight for PPR-Meta score |
| `abundance.primary_metric` | `tpm` | Primary abundance metric (`tpm`, `rpkm`, `mean`, `count`, `covered_fraction`) |
| `abundance.additional_metrics` | `[rpkm, mean, covered_fraction]` | Also computed, into `abundance_all_metrics.tsv` |
| `votu_clustering.enabled` | `true` | Dereplicate viral contigs into vOTUs and measure abundance on representatives |
| `votu_clustering.ani_threshold` | `95` | % identity for same-vOTU membership |
| `votu_clustering.target_coverage_threshold` | `85` | % of the absorbed contig the representative must cover |
| `test_mode` | `true` | Stub all tool calls with fast mocks (no databases needed) |
| `preprocessing.fastp.min_length` | `50` | fastp minimum read length |
| `preprocessing.fastp.qualified_quality` | `15` | fastp per-base quality threshold |
| `preprocessing.host_removal.enabled` | `false` | Remove host reads (minimap2); set `reference_genomes`. Recommended for tissue samples |
| `preprocessing.deduplication.enabled` | `false` | clumpify.sh read deduplication |
| `preprocessing.normalization.enabled` | `true` | bbnorm before assembly (assembly only; abundance uses un-normalized reads) |
| `assembly.mode` | `coassembly` | `coassembly` (pool all read samples, one assembly) or `per_sample` |
| `assembly.assembler` | `megahit` | `megahit` or `metaspades` |
| `assembly.min_contig_len` | `500` | Assembler minimum contig length |
| `assembly.quast` | `true` | QUAST report per assembly (included in MultiQC) |
| `virus_identification.tools.*` | see below | Enable/disable each detection tool |
| `virus_identification.n_chunks` | `1` | Split contigs into N chunks so geNomad/VirSorter2/VIBRANT/Cenote-Taker3 run as parallel jobs |
| `virus_identification.diamond.nr_scope` | `candidates` | Run DIAMOND vs nr on candidate contigs only (`candidates`) or all (`all`) |
| `virus_identification.diamond.bitscore_margin_pct` | `70` | Hits within this % of the best bitscore join the majority-vote call (placeholder — needs real-data tuning) |

**Detection tools** (`virus_identification.tools`): `genomad`, `virsorter2`, `vibrant`,
`cenote_taker3`, `diamond_refseq_viral` and `diamond_nr` are on by default;
`pprmeta` is off (phage/plasmid oriented, and it pulls in TensorFlow). VIBRANT
should stay on if phage filtering matters: the high-confidence phage rules key
on its calls.

Database paths are automatically set to `databases/` and require no editing.

---

## Databases

Databases are downloaded automatically the first time vine-mirus runs in
production mode (`test_mode: false`), into `databases/` at the repository root,
and reused by every later run.

| Database | Size | Tool | Download rule |
|----------|------|------|---------------|
| geNomad DB | ~3 GB | geNomad | `download_genomad_db` |
| VirSorter2 DB | ~10 GB | VirSorter2 | `download_virsorter2_db` |
| VIBRANT DB | ~11 GB | VIBRANT | `download_vibrant_db` |
| Cenote-Taker3 DBs | ~3 GB | Cenote-Taker3 | `download_cenote_taker3_db` |
| CheckV DB | ~3 GB | CheckV | `download_checkv_db` |
| NCBI taxdump | ~70 MB | taxopy + DIAMOND `makedb` | `download_taxdump` |
| `prot.accession2taxid` | ~2 GB | DIAMOND `makedb --taxonmap` | `download_prot_accession2taxid` |
| RefSeq viral DIAMOND | ~2 GB | DIAMOND | `build_diamond_refseq_viral_db` |
| nr DIAMOND | **~250 GB+** | DIAMOND | `build_diamond_nr_db` (opt-in) |
| PPR-Meta scripts | ~50 MB | PPR-Meta | `setup_pprmeta` |

Fetch everything the enabled tools need, without running any analysis:

```bash
snakemake --snakefile workflow/Snakefile databases \
          --config samples=config/samples.tsv test_mode=false \
          --cores 4 --use-conda
```

(The target name must come before `--config`; Snakemake otherwise reads it as
a config entry.)

### nr is opt-in

The nr search is Tier B's cellular-false-positive check, and nr is hundreds of
GB and many hours to download and build, so it is never fetched implicitly.
Pick one:

- `download_nr: true` in `config/config.yaml` to build it locally,
- point `diamond_nr_db` at an existing taxonomy-enabled `nr.dmnd` (e.g. a
  shared cluster path), or
- `virus_identification.tools.diamond_nr: false` to skip the search.

Until then, the `diamond_nr_db_missing` rule fails with exactly those three
options rather than silently starting a 250 GB download.

### Why each download is verified

Several of these downloaders exit 0 having written nothing usable, or having
written it one level deeper than the tool that consumes it then looks:

| Downloader | What it actually does |
|------------|----------------------|
| `genomad download-database X` | creates `X/genomad_db/` — `config.genomad_db` points at the nested path |
| `checkv download_database X` | creates `X/checkv-db-v1.N/` — resolved at run time, so no version is pinned here |
| `download-db.sh X` | creates `X/databases/`; VIBRANT also needs a sibling `files/`, which ships in the package. Passing the parent gives "could not identify KEGG HMM files" |
| `virsorter setup` | writes to the **current directory**, ignoring `--db-dir` and `--config`, while still reporting success. The rule cd's into the target first |

Every rule therefore checks for a specific expected file before its `.done`
sentinel is written, and `build_diamond_refseq_viral_db` takes RefSeq viral
protein parts until one 404s — fetching only `viral.1` silently truncates the
reference, and the searches then just find less.

The DIAMOND databases are built with `--taxonmap/--taxonnodes/--taxonnames`, so
every hit carries a real taxid. Joining accessions to taxids after the search
instead is what produced `staxids=0` on every hit in the original pipeline: a
well-formed result in which nothing could ever be classified as viral. The
DIAMOND rule fails if it ever sees that again.

### Hardware notes (CPUs without AVX2)

The reference deployment (`goldberg-r710`, 24 cores) has SSE4.2 but no AVX2, and
several current conda builds crash there with "Illegal instruction":

| Tool | Handling |
|------|----------|
| MEGAHIT | pinned to `1.2.8` in `envs/qc_assembly.yaml` (harmless on newer CPUs) |
| DIAMOND | `envs/checkv.yaml` pins `2.1.8`; the search envs use a current release |
| MMseqs2 | opt-in binary swap, see below |

MMseqs2 (inside Cenote-Taker3) has no suitable conda build, so
`envs/cenote_taker3.post-deploy.sh` downloads the static SSE4.1 build and swaps
it in. Snakemake runs post-deploy scripts whenever it creates the environment,
which is what keeps the fix from being silently reverted by the next env
rebuild. It is opt-in per machine:

```bash
export VINE_LEGACY_CPU=1   # before the run, on a CPU without AVX2
```

The script verifies the downloaded binary actually runs before replacing
anything, and keeps the conda build alongside it as `mmseqs.conda`.

---

## Running the Pipeline

### Local execution

A local profile is provided for the 24-core workstation. Unlike the Slurm
profile, its `cores` and `resources` are a budget shared by everything running
at once, so Snakemake packs jobs up to those limits — set `mem_mb` there to
roughly 90% of the machine's RAM.

```bash
# Production (downloads databases on first run)
snakemake --snakefile workflow/Snakefile \
          --config samples=config/my_samples.tsv test_mode=false \
          --profile profiles/local

# Test mode (mocks all tools — no databases, finishes in minutes)
snakemake --snakefile workflow/Snakefile \
          --config samples=config/my_samples.tsv \
          --cores 4 --use-conda --scheduler greedy
```

Boolean config values given on the command line are parsed explicitly, so
`test_mode=false`, `test_mode=False` and `test_mode=no` all work and anything
ambiguous is rejected. Snakemake itself evaluates `--config` values as Python,
where bare `false` is not a keyword and arrives as a non-empty string — which
would otherwise leave a "production" run quietly using the mocks.

Nested keys cannot be overridden with `--config` at all (`votu_clustering` is a
mapping, not a scalar). Copy `config/config.yaml`, edit it, and pass
`--configfile` instead.

### SLURM execution

A SLURM profile is provided in `profiles/slurm/`. Per-rule resource requirements
are defined in each rule's `resources:` block and passed automatically to the
SLURM executor.

```bash
snakemake --snakefile workflow/Snakefile \
          --config samples=config/my_samples.tsv test_mode=false \
          --profile profiles/slurm
```

Edit `profiles/slurm/config.yaml` to set your partition, account, and default
resource limits.

### Dry run (check workflow without executing)

```bash
snakemake --snakefile workflow/Snakefile \
          --config samples=config/my_samples.tsv -n
```

### Force re-run of specific rules

```bash
# Re-run everything from evidence scoring onward
snakemake --snakefile workflow/Snakefile \
          --config samples=config/my_samples.tsv \
          --forcerun evidence_score --cores 16 --use-conda
```

---

## Outputs

### Terminal outputs (results root)

| File | Description |
|------|-------------|
| `results/evidence_table.tsv` | **Reviewer table**: every contig screened, with each tool's call/score, the weighted `viral_score`, confidence tier, taxonomy, host domain, CheckV quality, and empty `reviewer_decision` / `reviewer_notes` columns |
| `results/phage_evidence_table.tsv` | High-confidence phage, same columns, for independent review |
| `results/viral_contigs.fa` | Candidate viral contigs (phage excluded unless `include_phage`) |
| `results/viral_contigs_abundance.tsv` | Per-contig summary with vOTU membership and one abundance column per sample |
| `results/qc/multiqc_report.html` | fastp + QUAST report |

**Evidence columns worth knowing:**

| Column | Meaning |
|--------|---------|
| `n_tools_called` / `n_tools_evaluated` | How many tools called this contig, out of how many ran |
| `viral_score` | Weighted mean of the tools that report a probability (geNomad, VirSorter2, VIBRANT, PPR-Meta), renormalized over the enabled ones. Cenote-Taker3 and the DIAMOND searches contribute calls, not scores |
| `confidence_tier` | 1/2/3 by tool agreement — **review priority only, never a filter** |
| `is_candidate` | Per `collation.candidate_rule` (`union` by default) |
| `taxonomy` / `taxonomy_source` | Lineage, and which tool it came from (geNomad → Cenote-Taker3 → DIAMOND) |
| `host_domain` | `prokaryote` / `eukaryote` / `unknown`, from lineage-clade lookup; unknown is retained |
| `cellular_hit_fraction` | Share of qualifying nr hits that look cellular. High values on a called contig flag the classic cellular false positive — reported, never auto-filtered |
| `high_confidence_phage`, `phage_rule_triggered` | Rule A (VIBRANT + prokaryote taxonomy) or Rule B (VIBRANT the only caller) |

### Two-tier design

Tier A runs the detection tools on every contig. Tier B (DIAMOND vs nr, CheckV)
runs only on the candidates Tier A produced, which is where most of the runtime
saving comes from. High-confidence phage is excluded from the candidate FASTA
before Tier B (unless `include_phage`), so the expensive searches skip it.

### Sanity checks

After parsing, the evidence step fails when an enabled tool matched **zero**
contig IDs, or when a call/score column is entirely constant. Both are silent
failures that return exit code 0 — three of six parsers in the original
pipeline matched nothing for exactly this reason. Set
`collation.strict_checks: false` to downgrade them to warnings.

### vOTU clustering and abundance

Viral contigs are dereplicated into **vOTUs** (viral operational taxonomic
units) at 95% ANI over 85% of the shorter sequence — the community convention,
and the method CheckV's `anicalc.py` / `aniclust.py` implement. Those scripts
live in the CheckV *repository*, not in the bioconda `checkv` package, so
`workflow/scripts/votu_cluster.py` implements the same method rather than
shelling out to something the conda env does not contain.

Clustering runs on `results/viral_contigs.fa`, i.e. **after** phage exclusion,
so no cluster can be orphaned by its representative being filtered out later.
CoverM then maps each sample's QC'd reads (never the bbnorm-normalized ones)
against the **representatives only**, with `--mapper minimap2-sr`. This keeps
one virus population's reads from being split across near-identical contigs,
and shrinks the mapping reference.

Abundance is therefore **per vOTU**: every member row in the summary reports
its cluster's abundance, and `is_votu_representative` marks the contig that
carried the mapping. Set `votu_clustering.enabled: false` for strictly
per-contig abundance.

`logs/votu_cluster.log` reports the singleton fraction and how much the set
collapsed. On a good co-assembly most redundancy is already resolved, so a high
singleton fraction is expected rather than a problem — but zero between-contig
alignments is also what an ID mismatch looks like, so the log says so
explicitly when that happens.

### Intermediate outputs

```
results/
├── samples/{sample|coasm}/
│   ├── qc/                  ← fastp output (also used for abundance mapping)
│   ├── assembly/            ← assembler output + QUAST
│   ├── contigs/             ← length-gated, canonically renamed contigs
│   └── <tool>/              ← per-tool output (+ chunks/ when scattered)
└── intermediate/
    ├── evidence/            ← evidence_tierA.tsv, candidates.fa
    ├── checkv/              ← quality_summary.tsv
    ├── votu/                ← votu_clusters.tsv, votu_representatives.fa
    └── abundance/           ← abundance_matrix.tsv (primary metric),
                               abundance_all_metrics.tsv (every metric)
```

## Test Mode

`test_mode: true` (default) replaces all tool calls with fast Python mocks that
generate random-but-realistic outputs. No databases are required, and the full
pipeline completes in seconds on any laptop.

Use test mode to:
- Verify the pipeline logic end-to-end before acquiring databases
- Develop and debug changes to the workflow
- Run CI/CD checks

```bash
# Explicit override to test mode
snakemake --snakefile workflow/Snakefile \
          --config samples=config/samples.tsv test_mode=true \
          --cores 4 --use-conda --scheduler greedy
```

Set `test_mode: false` in `config/config.yaml` (or pass `--config test_mode=false`)
for a real production run.

---

## Conda Environments

One environment per dependency island. There are more than the ideal four or
five because of genuine conflicts, noted below.

| Environment | Tools |
|-------------|-------|
| `envs/qc_assembly.yaml` | fastp, MEGAHIT (`=1.2.8`), metaSPAdes, bbmap, minimap2, samtools, seqkit |
| `envs/reports.yaml` | MultiQC, QUAST |
| `envs/genomad.yaml` | geNomad |
| `envs/virsorter2.yaml` | VirSorter2 (bioconda: `virsorter`) |
| `envs/vibrant.yaml` | VIBRANT `=1.2.1` (needs `numpy<1.24`) |
| `envs/cenote_taker3.yaml` | Cenote-Taker3 `=3.4.4` (+ an opt-in MMseqs2 swap, see Hardware notes) |
| `envs/diamond.yaml` | current DIAMOND, for the RefSeq viral and nr searches |
| `envs/checkv.yaml` | CheckV, BLAST (vOTU clustering), DIAMOND `=2.1.8` |
| `envs/pprmeta.yaml` | TensorFlow, biopython (PPR-Meta runtime) |
| `envs/viral_analysis.yaml` | CoverM, BWA, taxopy |
| `envs/scripts.yaml` | Python, pandas, biopython, taxopy (all pipeline scripts and mocks) |
| `envs/download.yaml` | curl, tar, unzip — the database download rules only |

`checkv` and `diamond` are separate because CheckV caps DIAMOND at `≤2.1.8`
while the nr and RefSeq searches want a current release. In test mode every
tool rule uses `envs/scripts.yaml` instead, because the mocks are Python and the
real tools are never invoked.

Environments are built automatically on first use by Snakemake (`--use-conda`).

---

## Citation

If you use vine-mirus in your research, please cite the individual tools:

- **geNomad**: Camargo *et al.*, Nature Biotechnology (2023)
- **VirSorter2**: Guo *et al.*, Microbiome (2021)
- **PPR-Meta**: Fang *et al.*, GigaScience (2019)
- **CheckV**: Nayfach *et al.*, Nature Biotechnology (2021)
- **DIAMOND**: Buchfink *et al.*, Nature Methods (2021)
- **CoverM**: github.com/wwood/CoverM
- **MEGAHIT**: Li *et al.*, Bioinformatics (2015)
- **fastp**: Chen *et al.*, Bioinformatics (2018)
