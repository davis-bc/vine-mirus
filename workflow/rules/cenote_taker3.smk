# Cenote-Taker3 — per chunk, then merged.
#
# Confirmed behaviours this rule works around (all from Jordan's live runs):
#   - -r/--run_title is a NAME (<18 chars, alnum/underscore), not a path;
#     Cenote-Taker3 creates {working_directory}/{run_title}/ itself.
#   - -db is a category selector (virion/rdrp/dnarep), NOT a database path;
#     the database goes to --cenote-dbs.
#   - it can exit 0 after printing "Exiting." and doing no work, so the real
#     summary file and its own completion marker are both checked.
#   - the params key must not be called "workdir": Snakemake reserves that.

_CT3_CMD = (
    "python {params.scripts}/mock_cenote_taker3.py "
    "{input.contigs} {output.summary} {params.min_len} 2> {log}"
    if TEST_MODE else
    """
    set -euo pipefail
    rm -rf {params.ct3_workdir}/{params.run_title}
    mkdir -p {params.ct3_workdir}
    cenotetaker3 -c {input.contigs} -r {params.run_title} \
        -p {params.prune} -wd {params.ct3_workdir} \
        -db {params.domains} --cenote-dbs {params.db} \
        -t {threads} &> {log}

    # Exit 0 is not proof of success for this tool: require its own completion
    # marker in the log as well as a real, non-empty summary file.
    if ! grep -q "Cenote-Taker finishing now" {log}; then
        echo "ERROR: cenotetaker3 exited without its completion marker" >> {log}
        exit 1
    fi
    s={params.ct3_workdir}/{params.run_title}/cenote_taker3_virus_summary.tsv
    if [ ! -s "$s" ]; then
        echo "ERROR: cenotetaker3 produced no virus summary at $s" >> {log}
        exit 1
    fi
    cp "$s" {output.summary}
    """
)

rule cenote_taker3:
    input:
        contigs="results/samples/{sample}/chunks/chunk_{chunk}.fa",
        db=[] if TEST_MODE else ["databases/cenote_taker3_db/.done"],
    output:
        summary="results/samples/{sample}/cenote_taker3/chunks/{chunk}.tsv",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=180,
    params:
        ct3_workdir="results/samples/{sample}/cenote_taker3/chunks",
        run_title="ct3_{chunk}",   # must stay <18 chars, [A-Za-z0-9_] only
        db=config["cenote_taker3_db"],
        domains=config["virus_identification"]["cenote_taker3"]["domains"],
        prune="True" if config["virus_identification"]["cenote_taker3"]["prune_prophage"] else "False",
        min_len=config["min_contig_length"],
        scripts=SCRIPTS,
    conda:
        # mocks need python/biopython, not the real tool
        "../../envs/scripts.yaml" if TEST_MODE else "../../envs/cenote_taker3.yaml"
    benchmark:
        "benchmarks/cenote_taker3/{sample}_{chunk}.tsv"
    log:
        "logs/cenote_taker3/{sample}_{chunk}.log"
    shell:
        _CT3_CMD


rule cenote_taker3_merge:
    input:
        chunks=expand("results/samples/{{sample}}/cenote_taker3/chunks/{chunk}.tsv", chunk=CHUNKS),
    output:
        summary="results/samples/{sample}/cenote_taker3/virus_summary.tsv",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15,
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/cenote_taker3/{sample}_merge.log"
    script:
        "../scripts/merge_tables.py"
