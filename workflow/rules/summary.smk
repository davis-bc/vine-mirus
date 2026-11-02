# Final per-contig summary: evidence table + vOTU membership + abundance, one
# row per viral contig, one abundance column per sample.

def _summary_input():
    inputs = {
        "evidence": "results/evidence_table.tsv",
        "viral_fa": "results/viral_contigs.fa",
        "abundance": "results/intermediate/abundance/abundance_matrix.tsv",
    }
    if VOTU:
        inputs["clusters"] = "results/intermediate/votu/votu_clusters.tsv"
    return inputs


rule summary:
    input:
        **_summary_input(),
    output:
        tsv="results/viral_contigs_abundance.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        samples=READS_SAMPLES,
        metric=PRIMARY_METRIC,
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/summary.log"
    script:
        "../scripts/summarize.py"
