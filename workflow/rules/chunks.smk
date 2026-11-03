# chunks.smk — split each contig set into N pieces so the slow virus-ID tools
# (geNomad, VirSorter2, VIBRANT, Cenote-Taker3) run as parallel jobs.
#
# virus_identification.n_chunks: 1 = no scatter (one chunk == the whole set).
# Each tool writes one output per chunk; a merge step concatenates them back to
# the single canonical table the evidence step reads, so scatter width never
# changes any downstream path.

rule split_contigs:
    input:
        contigs=get_contigs_path("{sample}"),
    output:
        # temp(): chunk files are intermediates, and deleting them means a
        # changed n_chunks can never leave stale chunks behind for the tools
        # to read (which silently truncates the contig set).
        chunks=temp(expand("results/samples/{{sample}}/chunks/chunk_{chunk}.fa", chunk=CHUNKS)),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        # Declared so Snakemake re-splits when the scatter width changes.
        n_chunks=N_CHUNKS,
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/split_contigs/{sample}.log"
    script:
        "../scripts/split_fasta.py"
