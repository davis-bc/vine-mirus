# VirSorter2 — per chunk, then merged.
# seqname values carry a "||full" / "||partial" / "||lt2gene" suffix that the
# evidence parser strips (normalize_contig_id).

_VIRSORTER2_CMD = (
    "python {params.scripts}/mock_virsorter2.py "
    "{input.contigs} {output.scores} {params.min_len} 2> {log}"
    if TEST_MODE else
    """
    set -euo pipefail
    rm -rf {params.outdir} && mkdir -p {params.outdir}
    virsorter run -i {input.contigs} -w {params.outdir} \
        --include-groups {params.groups} \
        --min-score {params.min_score} --min-length {params.min_len} \
        -j {threads} --db-dir {params.db} all 2> {log}
    cp {params.outdir}/final-viral-score.tsv {output.scores}
    """
)

rule virsorter2:
    input:
        contigs="results/samples/{sample}/chunks/chunk_{chunk}.fa",
        db=[] if TEST_MODE else ["databases/virsorter2_db/.done"],
    output:
        scores="results/samples/{sample}/virsorter2/chunks/{chunk}.tsv",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=180,
    params:
        outdir="results/samples/{sample}/virsorter2/chunks/{chunk}_run",
        db=config["virsorter2_db"],
        groups=",".join(config["virsorter2_groups"]),
        min_score=config["virus_identification"]["virsorter2_min_score"],
        min_len=config["min_contig_length"],
        scripts=SCRIPTS,
    conda:
        # mocks need python/biopython, not the real tool
        "../../envs/scripts.yaml" if TEST_MODE else "../../envs/virsorter2.yaml"
    benchmark:
        "benchmarks/virsorter2/{sample}_{chunk}.tsv"
    log:
        "logs/virsorter2/{sample}_{chunk}.log"
    shell:
        _VIRSORTER2_CMD


rule virsorter2_merge:
    input:
        chunks=expand("results/samples/{{sample}}/virsorter2/chunks/{chunk}.tsv", chunk=CHUNKS),
    output:
        scores="results/samples/{sample}/virsorter2/final-viral-score.tsv",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15,
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/virsorter2/{sample}_merge.log"
    script:
        "../scripts/merge_tables.py"
