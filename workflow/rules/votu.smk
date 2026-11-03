import os

# vOTU clustering — dereplicate the viral contigs into viral operational
# taxonomic units (95% ANI over 85% of the shorter sequence by default) and
# keep one representative per cluster.
#
# Clustering runs on results/viral_contigs.fa, i.e. AFTER phage exclusion, so
# a cluster can never be orphaned by its representative being filtered out
# later. Abundance is then measured against the representatives only, which is
# both cheaper and the standard way to avoid splitting one virus population's
# reads across near-identical contigs.
#
# Rules are only defined when votu_clustering.enabled; coverm.smk falls back to
# viral_contigs.fa as its reference when they are absent.

if VOTU:

    _BLASTN_CMD = (
        "python {params.scripts}/mock_blastn.py {input.viral_fa} {output.hits} 2> {log}"
        if TEST_MODE else
        """
        set -euo pipefail
        mkdir -p $(dirname {output.hits}) {params.dbdir}
        makeblastdb -in {input.viral_fa} -dbtype nucl -out {params.db} > {log} 2>&1
        blastn -task {params.task} \
               -query {input.viral_fa} -db {params.db} \
               -out {output.hits} \
               -outfmt "6 qseqid sseqid pident length qstart qend sstart send" \
               -evalue {params.evalue} \
               -max_target_seqs {params.max_target_seqs} \
               -num_threads {threads} >> {log} 2>&1
        """
    )

    rule votu_blastn:
        input:
            viral_fa="results/viral_contigs.fa",
            # Declared so a change to the mock re-runs the search: scripts
            # invoked through params.scripts are not otherwise tracked.
            lib=[f"{SCRIPTS}/mock_blastn.py"] if TEST_MODE else [],
        output:
            # temp(): a large all-vs-all table with no downstream use once the
            # clusters exist.
            hits=temp("results/intermediate/votu/all_vs_all.blastn.tsv"),
        threads: 8
        resources:
            mem_mb=16000,
            runtime=120,
        params:
            dbdir="results/intermediate/votu/blastdb",
            db="results/intermediate/votu/blastdb/viral_contigs",
            task=VOTU_CFG["blastn_task"],
            evalue=VOTU_CFG["blastn_evalue"],
            # Self-comparison of a whole contig set: the default of 500 would
            # silently truncate the neighbour list for abundant populations.
            max_target_seqs=VOTU_CFG["max_target_seqs"],
            # Derived from the tracked mock input when there is one; the real
            # search does not use params.scripts at all.
            scripts=lambda wc, input: os.path.dirname(input.lib[0]) if input.lib else SCRIPTS,
        conda:
            "../../envs/scripts.yaml" if TEST_MODE else "../../envs/checkv.yaml"
        benchmark:
            "benchmarks/votu_blastn.tsv"
        log:
            "logs/votu_blastn.log"
        shell:
            _BLASTN_CMD

    rule votu_cluster:
        input:
            hits="results/intermediate/votu/all_vs_all.blastn.tsv",
            viral_fa="results/viral_contigs.fa",
            # Tracked so edits to the clustering implementation re-run it.
            lib=f"{SCRIPTS}/votu_cluster.py",
        output:
            clusters="results/intermediate/votu/votu_clusters.tsv",
            reps="results/intermediate/votu/votu_representatives.fa",
        threads: 1
        resources:
            mem_mb=8000,
            runtime=60,
        params:
            min_ani=VOTU_CFG["ani_threshold"],
            min_tcov=VOTU_CFG["target_coverage_threshold"],
            min_qcov=VOTU_CFG["query_coverage_threshold"],
            # Derived from the tracked input so there is one source of truth.
            scripts=lambda wc, input: os.path.dirname(input.lib),
        conda:
            "../../envs/scripts.yaml"
        log:
            "logs/votu_cluster.log"
        shell:
            "python {params.scripts}/votu_cluster.py "
            "--blast {input.hits} --contigs {input.viral_fa} "
            "--min-ani {params.min_ani} --min-tcov {params.min_tcov} "
            "--min-qcov {params.min_qcov} "
            "--out-clusters {output.clusters} --out-reps {output.reps} "
            "> {log} 2>&1"
