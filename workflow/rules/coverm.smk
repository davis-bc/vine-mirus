# Abundance of the viral contigs across samples.
#
# The reference is the vOTU representatives when clustering is enabled, so
# reads from one virus population are not split across near-identical contigs.
# Mapping uses the QC'd reads, never the bbnorm-normalized ones.

def abundance_reference(wildcards=None):
    if VOTU:
        return "results/intermediate/votu/votu_representatives.fa"
    return "results/viral_contigs.fa"


if READS_SAMPLES:
    rule coverm:
        input:
            reference=abundance_reference,
            r1=expand("results/samples/{sample}/qc/{sample}_R1.fastq.gz", sample=READS_SAMPLES),
            r2=expand("results/samples/{sample}/qc/{sample}_R2.fastq.gz", sample=READS_SAMPLES),
        output:
            abundance="results/intermediate/abundance/abundance_matrix.tsv",
            all_metrics="results/intermediate/abundance/abundance_all_metrics.tsv",
        threads: 8
        resources:
            mem_mb=32000,
            runtime=120,
        params:
            samples=READS_SAMPLES,
            primary_metric=PRIMARY_METRIC,
            metrics=ABUNDANCE_METRICS,
            test_mode=TEST_MODE,
        conda:
            "../../envs/scripts.yaml" if TEST_MODE else "../../envs/viral_analysis.yaml"
        benchmark:
            "benchmarks/coverm.tsv"
        log:
            "logs/coverm.log"
        script:
            "../scripts/run_coverm.py"
else:
    # Contigs-only sample sheet: no reads to map, so the summary joins against
    # an empty matrix rather than the rule being absent from the DAG.
    rule coverm:
        output:
            abundance="results/intermediate/abundance/abundance_matrix.tsv",
            all_metrics="results/intermediate/abundance/abundance_all_metrics.tsv",
        log:
            "logs/coverm.log"
        shell:
            "mkdir -p $(dirname {output.abundance}) && "
            "printf 'contig_id\\n' > {output.abundance} && "
            "printf 'Contig\\n' > {output.all_metrics} && "
            "echo '[no reads] skipped CoverM' > {log}"
