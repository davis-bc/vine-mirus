# DIAMOND blastx searches.
#
# refseq_viral replaces Jordan's NCBI `blastx` against the same database:
# DIAMOND is orders of magnitude faster at comparable sensitivity for this
# step, and its makedb takes --taxonmap/--taxonnodes directly, which removes
# the fragile accession2taxid join that silently produced staxids=0 for every
# hit in his run.
#
# nr is the cellular-false-positive / host-domain check and is by far the most
# expensive search, so diamond.nr_scope defaults to "candidates": it only sees
# contigs some Tier A tool already called. The candidates FASTA is produced by
# the evidence step (evidence_score.smk).

_DMND = config["virus_identification"]["diamond"]

# 15 columns, as the evidence parser expects (bitscore at 12, staxids at 13)
DIAMOND_OUTFMT = ("6 qseqid sseqid pident length mismatch gapopen qstart qend "
                  "sstart send evalue bitscore staxids sscinames sskingdoms")

DIAMOND_DB = {
    "refseq_viral": config["diamond_refseq_viral_db"],
    "nr": config["diamond_nr_db"],
}


def diamond_query(wildcards):
    """refseq_viral screens every contig; nr screens candidates by default."""
    if wildcards.db == "nr" and _DMND["nr_scope"] == "candidates":
        return "results/intermediate/evidence/candidates.fa"
    return get_contigs_path(wildcards.sample)


_DIAMOND_CMD = (
    "python {params.scripts}/mock_diamond.py "
    "{input.query} {output.hits} {params.min_len} 2> {log}"
    if TEST_MODE else
    """
    set -euo pipefail
    mkdir -p $(dirname {output.hits})
    diamond blastx -q {input.query} -d {params.db} -o {output.hits} \
        -e {params.evalue} --outfmt {params.outfmt} \
        --max-target-seqs {params.max_target} {params.sensitivity} \
        --threads {threads} 2> {log}

    # A well-formed file of hits that all carry staxids=0 means the database
    # was built without usable taxonomy: the majority-vote call would then
    # classify nothing as viral. Fail loudly instead.
    if [ -s {output.hits} ]; then
        nz=$(cut -f13 {output.hits} | grep -cvE '^(0|)$' || true)
        if [ "$nz" -eq 0 ]; then
            echo "ERROR: every hit has staxids=0 - {params.db} lacks taxonomy" >> {log}
            exit 1
        fi
    fi
    """
)

rule diamond:
    input:
        query=diamond_query,
        db=[] if TEST_MODE else lambda wc: [DIAMOND_DB[wc.db]],
    output:
        hits="results/samples/{sample}/diamond_{db}/hits.tsv",
    wildcard_constraints:
        db="refseq_viral|nr",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=240,
    params:
        db=lambda wc: DIAMOND_DB[wc.db],
        evalue=_DMND["evalue"],
        outfmt=DIAMOND_OUTFMT,
        max_target=_DMND["max_target_seqs"],
        sensitivity=_DMND["sensitivity"],
        min_len=config["min_contig_length"],
        scripts=SCRIPTS,
    conda:
        # mocks need python/biopython, not the real tool
        "../../envs/scripts.yaml" if TEST_MODE else "../../envs/diamond.yaml"
    benchmark:
        "benchmarks/diamond/{sample}_{db}.tsv"
    log:
        "logs/diamond/{sample}_{db}.log"
    shell:
        _DIAMOND_CMD
