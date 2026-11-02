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
       geNomad      VirSorter2    PPR-Meta
          └─────────────┼─────────────┘
                        ▼
               evidence_score  ←  weighted fusion
                        │
                        ▼
            taxonomy_checkv  ←  DIAMOND + NCBI taxdump + CheckV
                        │
                        ▼
                     CoverM  ←  read-level abundance per contig
                        │
                        ▼
          results/viral_contigs_abundance.tsv
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
| `coverm_metric` | `rpkm` | CoverM abundance metric (`rpkm`, `tpm`, `mean`) |
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

Database paths are automatically set to `databases/` and require no editing.

---

## Databases

Databases are downloaded automatically the first time vine-mirus runs in
production mode (`test_mode: false`). They are stored in `databases/` at the
repository root and reused across all subsequent runs.

| Database | Size | Tool | Download rule |
|----------|------|------|---------------|
| geNomad DB | ~3 GB | geNomad | `download_genomad_db` |
| VirSorter2 DB | ~10 GB | VirSorter2 | `download_virsorter2_db` |
| CheckV DB | ~3 GB | CheckV | `download_checkv_db` |
| NCBI taxdump | ~70 MB | DIAMOND/taxopy | `download_taxdump` |
| Viral RefSeq DIAMOND | ~2 GB | DIAMOND | `build_diamond_db` |
| PPR-Meta scripts | ~50 MB | PPR-Meta | `setup_pprmeta` |

To trigger all database downloads without running the full pipeline:

```bash
snakemake --snakefile workflow/Snakefile \
          --config samples=config/samples.tsv test_mode=false \
          --until download_genomad_db download_virsorter2_db \
                  download_checkv_db download_taxdump build_diamond_db \
          --cores 4 --use-conda
```

---

## Running the Pipeline

### Local execution

```bash
# Production (downloads databases on first run)
snakemake --snakefile workflow/Snakefile \
          --config samples=config/my_samples.tsv test_mode=false \
          --cores 16 --use-conda

# Test mode (mocks all tools — no databases, finishes in minutes)
snakemake --snakefile workflow/Snakefile \
          --config samples=config/my_samples.tsv \
          --cores 4 --use-conda --scheduler greedy
```

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
| `results/viral_contigs_abundance.tsv` | Main output: viral contigs with taxonomy, completeness, and per-sample RPKM |
| `results/viral_contigs.fa` | FASTA of all viral contigs passing the evidence threshold |

**`viral_contigs_abundance.tsv` columns:**

| Column | Description |
|--------|-------------|
| `contig_id` | Contig identifier from assembly |
| `sample_of_origin` | Sample in which the contig was assembled |
| `length` | Contig length (bp) |
| `viral_score` | Weighted evidence score (0–1) |
| `taxonomy` | NCBI lineage: `superkingdom;phylum;class;order;family;genus;species` |
| `checkv_quality` | CheckV completeness: Complete / High-quality / Medium-quality / Low-quality / Not-determined |
| `<sample> RPKM` | One column per sample: reads per kilobase per million mapped reads |

### Intermediate outputs

All intermediate files are nested under `results/` to keep the root uncluttered:

```
results/
├── viral_contigs_abundance.tsv       ← terminal
├── viral_contigs.fa                  ← terminal
├── samples/
│   └── {sample}/
│       ├── qc/                       ← fastp trimmed reads + reports
│       ├── assembly/                 ← MEGAHIT contigs
│       ├── genomad/                  ← geNomad virus summary
│       ├── virsorter2/               ← VirSorter2 scores
│       └── pprmeta/                  ← PPR-Meta predictions
└── intermediate/
    ├── evidence/viral_scores.tsv     ← fused per-contig scores
    ├── taxonomy/viral_taxonomy.tsv   ← DIAMOND taxonomy hits
    ├── checkv/quality_summary.tsv    ← CheckV completeness
    └── abundance/abundance_matrix.tsv ← CoverM raw matrix
```

---

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

vine-mirus uses 6 minimal conda environments:

| Environment | Tools |
|-------------|-------|
| `envs/qc_assembly.yaml` | fastp, MEGAHIT |
| `envs/genomad.yaml` | geNomad |
| `envs/virsorter2.yaml` | VirSorter2 (bioconda: `virsorter`) |
| `envs/pprmeta.yaml` | TensorFlow, biopython (PPR-Meta runtime) |
| `envs/viral_analysis.yaml` | CheckV, DIAMOND, CoverM, BWA, taxopy |
| `envs/scripts.yaml` | Python, pandas, biopython (evidence scoring + summary) |

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
