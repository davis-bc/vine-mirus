# VIBRANT — per chunk, then merged.
#
# Two things Jordan confirmed the hard way:
#   - VIBRANT_run.py needs TWO database paths: -d <vibrant_db>/databases and
#     -m its sibling <vibrant_db>/files. Passing only -d fails with
#     "could not identify KEGG HMM files".
#   - the genome-quality table's "scaffold" column holds the FULL original
#     FASTA header, and its quality column is "Quality" (capital Q).

_VIBRANT_CMD = (
    "python {params.scripts}/mock_vibrant.py "
    "{input.contigs} {output.quality} {params.min_len} 2> {log}"
    if TEST_MODE else
    """
    set -euo pipefail
    rm -rf {params.outdir} && mkdir -p {params.outdir}
    base=$(basename {input.contigs} .fa)
    VIBRANT_run.py -i {input.contigs} -folder {params.outdir} \
        -d {params.db}/databases -m {params.db}/files \
        -t {threads} 2> {log}
    q={params.outdir}/VIBRANT_${{base}}/VIBRANT_results_${{base}}/VIBRANT_genome_quality_${{base}}.tsv
    if [ ! -s "$q" ]; then
        echo "ERROR: VIBRANT produced no genome-quality table at $q" >> {log}
        exit 1
    fi
    cp "$q" {output.quality}
    """
)

rule vibrant:
    input:
        contigs="results/samples/{sample}/chunks/chunk_{chunk}.fa",
        db=[] if TEST_MODE else ["databases/vibrant_db/.done"],
    output:
        quality="results/samples/{sample}/vibrant/chunks/{chunk}.tsv",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=180,
    params:
        outdir="results/samples/{sample}/vibrant/chunks/{chunk}_run",
        db=config["vibrant_db"],
        min_len=config["min_contig_length"],
        scripts=SCRIPTS,
    conda:
        # mocks need python/biopython, not the real tool
        "../../envs/scripts.yaml" if TEST_MODE else "../../envs/vibrant.yaml"
    benchmark:
        "benchmarks/vibrant/{sample}_{chunk}.tsv"
    log:
        "logs/vibrant/{sample}_{chunk}.log"
    shell:
        _VIBRANT_CMD


rule vibrant_merge:
    input:
        chunks=expand("results/samples/{{sample}}/vibrant/chunks/{chunk}.tsv", chunk=CHUNKS),
    output:
        quality="results/samples/{sample}/vibrant/genome_quality.tsv",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15,
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/vibrant/{sample}_merge.log"
    script:
        "../scripts/merge_tables.py"
