# PPR-Meta (optional; virus_identification.tools.pprmeta).
# Disabled by default: phage/plasmid oriented, and it pulls in TensorFlow.

PPRMETA_DIR = "databases/pprmeta"
PPRMETA_SCRIPT = f"{PPRMETA_DIR}/PPR_Meta.py"

_PPRMETA_CMD = (
    "python {params.scripts}/mock_pprmeta.py "
    "{input.contigs} {output.results} {params.min_len} 2> {log}"
    if TEST_MODE else
    """
    set -euo pipefail
    mkdir -p {params.outdir}
    python {params.scripts}/filter_fasta_by_length.py \
        {input.contigs} {params.outdir}/filtered.fa {params.min_len}
    python3 {input.script} {params.outdir}/filtered.fa {output.results} 2> {log}
    """
)


rule setup_pprmeta:
    r"""Download PPR-Meta v1.1 from GitHub into databases/pprmeta/."""
    output:
        script=PPRMETA_SCRIPT,
    params:
        url="https://github.com/zhenchengfang/PPR-Meta/releases/download/v1.1/PPR_Meta_v_1_1.zip",
        outdir=PPRMETA_DIR,
    log:
        "logs/setup_pprmeta.log"
    shell:
        r"""
        mkdir -p {params.outdir}
        curl -L {params.url} -o {params.outdir}/PPR_Meta.zip 2> {log}
        unzip -o {params.outdir}/PPR_Meta.zip -d {params.outdir} >> {log} 2>&1
        find {params.outdir} -name "PPR_Meta.py" ! -path "{output.script}" \
            -exec mv {{}} {output.script} \; 2>> {log} || true
        """


rule pprmeta:
    input:
        contigs=get_contigs,
        **({} if TEST_MODE else {"script": PPRMETA_SCRIPT}),
    output:
        results="results/samples/{sample}/pprmeta/pprmeta.csv",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=60,
    params:
        outdir="results/samples/{sample}/pprmeta",
        min_len=config["min_contig_length"],
        scripts=SCRIPTS,
    conda:
        # mocks need python/biopython, not the real tool
        "../../envs/scripts.yaml" if TEST_MODE else "../../envs/pprmeta.yaml"
    benchmark:
        "benchmarks/pprmeta/{sample}.tsv"
    log:
        "logs/pprmeta/{sample}.log"
    shell:
        _PPRMETA_CMD
