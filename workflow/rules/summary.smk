# Final per-contig summary: evidence table + abundance, one row per viral
# contig, one abundance column per sample.

rule summary:
    input:
        evidence="results/evidence_table.tsv",
        viral_fa="results/viral_contigs.fa",
        abundance="results/intermediate/abundance/abundance_matrix.tsv",
    output:
        tsv="results/viral_contigs_abundance.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        samples=READS_SAMPLES,
        metric=config["coverm_metric"],
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/summary.log"
    script:
        "../scripts/summarize.py"
