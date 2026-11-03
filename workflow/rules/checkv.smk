# CheckV completeness/contamination — Tier B, so it runs on candidate contigs
# only rather than every assembled contig.
#
# The database downloader creates a nested versioned directory, so the real
# path is databases/checkv_db/checkv-db-v1.x (resolved in databases.smk).

_CHECKV_CMD = (
    "python {params.scripts}/mock_checkv.py "
    "{input.candidates} {output.summary} 2> {log}"
    if TEST_MODE else
    """
    set -euo pipefail
    rm -rf {params.outdir} && mkdir -p {params.outdir}
    db=$(find {params.db} -maxdepth 1 -name 'checkv-db-*' -type d | head -1)
    if [ -z "$db" ]; then db={params.db}; fi
    checkv end_to_end {input.candidates} {params.outdir} -d "$db" -t {threads} 2> {log}
    cp {params.outdir}/quality_summary.tsv {output.summary}
    """
)

rule checkv:
    input:
        candidates="results/intermediate/evidence/candidates.fa",
        db=[] if TEST_MODE else ["databases/checkv_db/.done"],
    output:
        summary="results/intermediate/checkv/quality_summary.tsv",
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120,
    params:
        outdir="results/intermediate/checkv/run",
        db=config["checkv_db"],
        scripts=SCRIPTS,
    conda:
        "../../envs/scripts.yaml" if TEST_MODE else "../../envs/checkv.yaml"
    benchmark:
        "benchmarks/checkv.tsv"
    log:
        "logs/checkv.log"
    shell:
        _CHECKV_CMD
