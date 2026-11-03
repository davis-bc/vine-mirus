# geNomad — per chunk, then merged to one virus_summary.tsv per contig set.

_GENOMAD_CMD = (
    "python {params.scripts}/mock_genomad.py "
    "{input.contigs} {output.summary} {params.min_len} 2> {log}"
    if TEST_MODE else
    """
    set -euo pipefail
    rm -rf {params.outdir} && mkdir -p {params.outdir}
    genomad end-to-end --min-score 0.0 --min-length {params.min_len} \
        --threads {threads} --cleanup \
        {input.contigs} {params.outdir} {params.db} 2> {log}
    # geNomad names its output after the input file's basename
    base=$(basename {input.contigs} .fa)
    cp {params.outdir}/${{base}}_summary/${{base}}_virus_summary.tsv {output.summary}
    """
)

rule genomad:
    input:
        contigs="results/samples/{sample}/chunks/chunk_{chunk}.fa",
        db=[] if TEST_MODE else ["databases/genomad_db/.done"],
    output:
        summary="results/samples/{sample}/genomad/chunks/{chunk}.tsv",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        outdir="results/samples/{sample}/genomad/chunks/{chunk}_run",
        db=config["genomad_db"],
        min_len=config["min_contig_length"],
        scripts=SCRIPTS,
    conda:
        # mocks need python/biopython, not the real tool
        "../../envs/scripts.yaml" if TEST_MODE else "../../envs/genomad.yaml"
    benchmark:
        "benchmarks/genomad/{sample}_{chunk}.tsv"
    log:
        "logs/genomad/{sample}_{chunk}.log"
    shell:
        _GENOMAD_CMD


rule genomad_merge:
    input:
        chunks=expand("results/samples/{{sample}}/genomad/chunks/{chunk}.tsv", chunk=CHUNKS),
    output:
        summary="results/samples/{sample}/genomad/virus_summary.tsv",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15,
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/genomad/{sample}_merge.log"
    script:
        "../scripts/merge_tables.py"
