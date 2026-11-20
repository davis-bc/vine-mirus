"""Every setting `vine-mirus config` shows: where it lives, its type, and help.

One entry per key. A new key in config/config.yaml (or the execution block in
project.py) needs an entry here too, or the TUI will not show it.

Kinds:
  bool, int, float      - plain scalars (int/float take an optional lo/hi range)
  str, path             - free text (path is shown with a path hint)
  evalue                - a number kept as a STRING in YAML ("1e-5"), as the rules expect
  choice                - one of `choices`
  list                  - comma-separated; `choices` (if given) restricts the items
  words                 - space-separated string; `choices` restricts the words
"""

from __future__ import annotations

from dataclasses import dataclass, field

# (key, glyph, title) - the side-panel tabs, in order. Number keys 1-7 switch.
TABS = [
    ("run", "≋", "Run"),
    ("samples", "><>", "Samples & reads"),
    ("assembly", "▲", "Assembly"),
    ("detection", "⚓", "Detection"),
    ("evidence", "⚖", "Evidence"),
    ("votu", "◎", "vOTUs & abundance"),
    ("databases", "⛁", "Databases"),
]

TOOLS = ("genomad", "virsorter2", "vibrant", "cenote_taker3", "diamond_refseq_viral", "diamond_nr", "pprmeta")


@dataclass(frozen=True)
class Field:
    tab: str
    group: str
    path: tuple[str, ...]
    label: str
    kind: str
    help: str
    choices: tuple[str, ...] = ()
    lo: float | None = None
    hi: float | None = None
    show: str | None = None          # "local" / "slurm": only shown for that executor
    extra: dict = field(default_factory=dict)

    @property
    def widget_id(self) -> str:
        return "f-" + "--".join(self.path)

    @property
    def dotted(self) -> str:
        return ".".join(self.path)


def _f(tab, group, path, label, kind, help, **kw) -> Field:
    return Field(tab, group, tuple(path.split(".")), label, kind, help, **kw)


FIELDS: list[Field] = [
    # ── Run ───────────────────────────────────────────────────────────────────
    _f("run", "Where it runs", "execution.executor", "Executor", "choice",
       "local runs jobs on this machine; slurm submits each job to the cluster scheduler.",
       choices=("local", "slurm")),
    _f("run", "Where it runs", "test_mode", "Test mode", "bool",
       "Mock every tool (no databases, seconds). Turn OFF for a real run - left on, a 'production' run silently uses mocks."),
    _f("run", "Local budget", "execution.cores", "Total cores", "int",
       "Cores shared by all running jobs. Snakemake packs jobs up to this limit.", lo=1, show="local"),
    _f("run", "Local budget", "execution.mem_mb", "Total memory (MB)", "int",
       "Memory budget shared by all running jobs; set ~90% of RAM. Heavy rules ask 32000 each.", lo=1000, show="local"),
    _f("run", "Slurm", "execution.slurm_account", "Account", "str",
       "Account to charge (--account). Blank uses the cluster default.", show="slurm"),
    _f("run", "Slurm", "execution.slurm_partition", "Default partition", "str",
       "Partition for every rule that does not request its own. Check `sinfo` for names.", show="slurm"),
    _f("run", "Slurm", "execution.slurm_mem_mb", "Default memory (MB)", "int",
       "Memory for rules that declare none; the heavy rules set their own.", lo=500, show="slurm"),
    _f("run", "Slurm", "execution.slurm_runtime", "Default walltime (min)", "int",
       "Walltime for rules that declare none.", lo=1, show="slurm"),
    _f("run", "Slurm", "execution.jobs", "Max queued jobs", "int",
       "Most jobs submitted at once.", lo=1, show="slurm"),
    _f("run", "Behaviour", "execution.keep_going", "Keep going after a failure", "bool",
       "Finish independent jobs when one fails, instead of stopping the whole run."),
    _f("run", "Behaviour", "execution.rerun_incomplete", "Re-run incomplete jobs", "bool",
       "Redo jobs an interrupted run left half-written."),
    _f("run", "Behaviour", "execution.latency_wait", "Filesystem latency (s)", "int",
       "Seconds to wait for outputs to appear on a shared/network filesystem.", lo=0),
    _f("run", "Shared stores", "execution.conda_prefix", "Conda env store", "path",
       "Where tool envs are built, shared by every project so they are solved once."),
    _f("run", "Shared stores", "execution.databases_dir", "Database store", "path",
       "Shared reference databases (tens to hundreds of GB). <project>/databases links here."),
    _f("run", "Shared stores", "allow_no_conda", "Allow running without conda", "bool",
       "Only for setups with every tool already on PATH (e.g. a container). vine-mirus always enables conda."),

    # ── Samples & reads ───────────────────────────────────────────────────────
    _f("samples", "Sample sheet", "samples", "Sample sheet", "path",
       "TSV with a 'sample' column plus r1/r2 (reads) and/or contigs. Relative paths resolve from this project directory."),
    _f("samples", "fastp", "preprocessing.fastp.min_length", "Min read length (bp)", "int",
       "Reads shorter than this after trimming are dropped.", lo=1),
    _f("samples", "fastp", "preprocessing.fastp.qualified_quality", "Qualified base quality", "int",
       "Phred score a base needs to count as qualified (Jordan's pipeline used 30).", lo=0, hi=60),
    _f("samples", "Host removal", "preprocessing.host_removal.enabled", "Remove host reads", "bool",
       "minimap2 against the host genome(s); keeps pairs where both mates are unmapped. Strongly recommended for fish tissue."),
    _f("samples", "Host removal", "preprocessing.host_removal.reference_genomes", "Host reference FASTA(s)", "list",
       "Comma-separated paths; concatenated and indexed once. Required when host removal is on."),
    _f("samples", "Deduplication", "preprocessing.deduplication.enabled", "Deduplicate reads", "bool",
       "clumpify.sh dedupe=t before assembly and mapping."),
    _f("samples", "Normalization", "preprocessing.normalization.enabled", "Normalize before assembly", "bool",
       "bbnorm digital normalization. Affects assembly input only; abundance always uses un-normalized reads."),
    _f("samples", "Normalization", "preprocessing.normalization.target_depth", "Target depth", "int",
       "bbnorm target k-mer depth.", lo=1),
    _f("samples", "Normalization", "preprocessing.normalization.min_depth", "Min depth", "int",
       "bbnorm drops reads below this depth (Jordan: benchmark 0 vs 2 on real data).", lo=0),

    # ── Assembly ──────────────────────────────────────────────────────────────
    _f("assembly", "Assembly", "assembly.mode", "Mode", "choice",
       "coassembly pools all read samples into one assembly (contigs 'coasm__'); per_sample assembles each ('{sample}__').",
       choices=("coassembly", "per_sample")),
    _f("assembly", "Assembly", "assembly.assembler", "Assembler", "choice",
       "MEGAHIT is fast and memory-light; metaSPAdes is slower and hungrier.", choices=("megahit", "metaspades")),
    _f("assembly", "Assembly", "assembly.min_contig_len", "Assembler min contig (bp)", "int",
       "Assembler pre-filter; the viral tools gate separately on the length below.", lo=100),
    _f("assembly", "Assembly", "assembly.quast", "QUAST report", "bool",
       "QUAST on each assembly, included in the MultiQC report."),
    _f("assembly", "Contig gate", "min_contig_length", "Min contig for detection (bp)", "int",
       "Contigs shorter than this never reach the viral detection tools.", lo=100),

    # ── Detection ─────────────────────────────────────────────────────────────
    *[_f("detection", "Tools", f"virus_identification.tools.{t}", label, "bool", help)
      for t, label, help in [
          ("genomad", "geNomad", "Marker + neural-net classifier; first choice for taxonomy."),
          ("virsorter2", "VirSorter2", "Multi-classifier viral scoring across the groups below."),
          ("vibrant", "VIBRANT", "Phage-oriented. Keep on: the phage-exclusion rules depend on its calls."),
          ("cenote_taker3", "Cenote-Taker3", "Hallmark-gene detection, incl. RdRp - key for eukaryotic RNA viruses."),
          ("diamond_refseq_viral", "DIAMOND vs RefSeq viral", "blastx against viral proteins; adds a call and a taxonomy fallback."),
          ("diamond_nr", "DIAMOND vs nr (Tier B)", "Candidates only: flags cellular false positives. nr is ~250 GB."),
          ("pprmeta", "PPR-Meta", "Phage/plasmid classifier; off by default (poor fit for eukaryotic viruses, needs TensorFlow)."),
      ]],
    _f("detection", "Scatter", "virus_identification.n_chunks", "Contig chunks", "int",
       "Split contigs into this many parallel jobs for the slow Tier A tools. 1 = no splitting.", lo=1),
    _f("detection", "Call thresholds", "virus_identification.genomad_min_score", "geNomad min score", "float",
       "virus_score at or above this counts as a geNomad call.", lo=0, hi=1),
    _f("detection", "Call thresholds", "virus_identification.virsorter2_min_score", "VirSorter2 min score", "float",
       "max_score at or above this counts as a VirSorter2 call (also passed to VirSorter2).", lo=0, hi=1),
    _f("detection", "Call thresholds", "virus_identification.pprmeta_min_score", "PPR-Meta min score", "float",
       "Phage_score at or above this counts as a PPR-Meta call.", lo=0, hi=1),
    _f("detection", "VirSorter2", "virsorter2_groups", "Viral groups", "list",
       "Classifier groups to run. Keep ssDNA, NCLDV and RNA for eukaryotic viruses.",
       choices=("dsDNAphage", "ssDNA", "NCLDV", "RNA", "lavidaviridae")),
    _f("detection", "Cenote-Taker3", "virus_identification.cenote_taker3.domains", "Hallmark databases", "words",
       "Space-separated: virion, rdrp, dnarep. All three for broad eukaryotic-virus detection.",
       choices=("virion", "rdrp", "dnarep")),
    _f("detection", "Cenote-Taker3", "virus_identification.cenote_taker3.prune_prophage", "Prune prophages", "bool",
       "Trim host flanks from prophage calls."),
    _f("detection", "DIAMOND", "virus_identification.diamond.evalue", "E-value", "evalue",
       "Maximum e-value for a hit.", lo=0),
    _f("detection", "DIAMOND", "virus_identification.diamond.max_target_seqs", "Max target seqs", "int",
       "Hits kept per contig.", lo=1),
    _f("detection", "DIAMOND", "virus_identification.diamond.top_n_hits", "Top hits for the vote", "int",
       "Hits considered in the majority-vote host call.", lo=1),
    _f("detection", "DIAMOND", "virus_identification.diamond.bitscore_margin_pct", "Bitscore margin (%)", "float",
       "Hits within this % of the best bitscore join the vote. Placeholder pending real-data benchmarking.", lo=0, hi=100),
    _f("detection", "DIAMOND", "virus_identification.diamond.sensitivity", "Sensitivity", "choice",
       "DIAMOND sensitivity mode; slower modes find more distant homologs.",
       choices=("--fast", "--mid-sensitive", "--sensitive", "--more-sensitive", "--very-sensitive", "--ultra-sensitive")),
    _f("detection", "DIAMOND", "virus_identification.diamond.nr_scope", "nr scope", "choice",
       "candidates searches nr only for contigs a Tier A tool called (the big runtime saving); all searches every contig.",
       choices=("candidates", "all")),

    # ── Evidence ──────────────────────────────────────────────────────────────
    _f("evidence", "Candidate rule", "collation.candidate_rule", "Candidate rule", "choice",
       "union: any tool called it (high recall). score: weighted viral_score >= threshold. both: must satisfy both.",
       choices=("union", "score", "both")),
    _f("evidence", "Candidate rule", "viral_score_threshold", "viral_score threshold", "float",
       "Used by the score/both rules. viral_score is always reported either way.", lo=0, hi=1),
    *[_f("evidence", "Score weights", f"tool_weights.{t}", f"{name} weight", "float",
         "Weight in viral_score; renormalized over the enabled tools that report a score.", lo=0)
      for t, name in [("genomad", "geNomad"), ("virsorter2", "VirSorter2"), ("vibrant", "VIBRANT"), ("pprmeta", "PPR-Meta")]],
    _f("evidence", "Review tiers (never filter)", "collation.confidence_tiers.tier_1_min_tools", "Tier 1: min tools", "int",
       "Review priority only - tiers never remove a contig.", lo=1),
    _f("evidence", "Review tiers (never filter)", "collation.confidence_tiers.tier_2_min_tools", "Tier 2: min tools", "int",
       "Must be <= tier 1.", lo=1),
    _f("evidence", "Review tiers (never filter)", "collation.confidence_tiers.tier_3_min_tools", "Tier 3: min tools", "int",
       "Must be <= tier 2.", lo=1),
    _f("evidence", "Safety checks", "collation.strict_checks", "Strict checks", "bool",
       "Fail when an enabled tool matched zero contigs or a column is constant - the signature of a silent ID break. Keep on."),
    _f("evidence", "Phage", "phage_filtering.include_phage", "Include phage", "bool",
       "Off: high-confidence phage skips Tier B and the main outputs (it is always written to its own table)."),

    # ── vOTUs & abundance ─────────────────────────────────────────────────────
    _f("votu", "vOTU clustering", "votu_clustering.enabled", "Cluster into vOTUs", "bool",
       "Dereplicate viral contigs; abundance is measured on one representative per vOTU."),
    _f("votu", "vOTU clustering", "votu_clustering.ani_threshold", "ANI (%)", "float",
       "Community convention: 95% ANI.", lo=0, hi=100),
    _f("votu", "vOTU clustering", "votu_clustering.target_coverage_threshold", "Target coverage (%)", "float",
       "Share of the absorbed contig the representative must cover. Convention: 85%.", lo=0, hi=100),
    _f("votu", "vOTU clustering", "votu_clustering.query_coverage_threshold", "Query coverage (%)", "float",
       "0 = no requirement (aniclust default).", lo=0, hi=100),
    _f("votu", "vOTU clustering", "votu_clustering.blastn_evalue", "blastn e-value", "evalue",
       "Strict, since this is a same-dataset self-comparison.", lo=0),
    _f("votu", "vOTU clustering", "votu_clustering.blastn_task", "blastn task", "choice",
       "megablast suits highly similar sequences and is much faster.", choices=("megablast", "dc-megablast", "blastn")),
    _f("votu", "vOTU clustering", "votu_clustering.max_target_seqs", "Max target seqs", "int",
       "Neighbours kept per contig; blastn's default 500 would truncate lists.", lo=1),
    _f("votu", "Abundance", "abundance.primary_metric", "Primary metric", "choice",
       "Written to the summary. TPM corrects for length and depth, so it compares best across samples.",
       choices=("tpm", "rpkm", "mean", "count", "covered_fraction")),
    _f("votu", "Abundance", "abundance.additional_metrics", "Additional metrics", "list",
       "Comma-separated; written to abundance_all_metrics.tsv for review.",
       choices=("tpm", "rpkm", "mean", "count", "covered_fraction")),

    # ── Databases ─────────────────────────────────────────────────────────────
    _f("databases", "nr", "download_nr", "Download nr", "bool",
       "Fetch and build NCBI nr (~250 GB download, hours). Off: provide diamond_nr_db yourself or disable DIAMOND vs nr."),
    *[_f("databases", "Paths (relative to the project)", key, label, "path", help)
      for key, label, help in [
          ("genomad_db", "geNomad", "Points one level inside the download dir (genomad_db/genomad_db)."),
          ("virsorter2_db", "VirSorter2", "virsorter setup target."),
          ("checkv_db", "CheckV", "The versioned subdirectory is resolved at run time."),
          ("vibrant_db", "VIBRANT", "Parent of VIBRANT's databases/ and files/."),
          ("cenote_taker3_db", "Cenote-Taker3", "Passed as --cenote-dbs."),
          ("taxdump_dir", "NCBI taxdump", "nodes.dmp + names.dmp, for DIAMOND taxonomy."),
          ("diamond_refseq_viral_db", "DIAMOND RefSeq viral", "Built with taxonomy so hits carry taxids."),
          ("diamond_nr_db", "DIAMOND nr", "Set to an existing nr.dmnd to skip the download."),
      ]],
]

FIELD_BY_PATH = {f.path: f for f in FIELDS}
