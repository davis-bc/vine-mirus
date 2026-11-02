rule prep_contigs:
    """Length-gate contigs once (min_contig_length) and rename to '{sample}__<id>'.

    Every viral-ID tool and every join downstream uses these IDs, so contigs
    from different samples can no longer collide (MEGAHIT reuses k141_N).
    """
    input:
        contigs=get_raw_contigs,
    output:
        contigs="results/samples/{sample}/contigs/{sample}_contigs.fa",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        prefix="{sample}",
        min_len=config["min_contig_length"],
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/prep_contigs/{sample}.log"
    script:
        "../scripts/prep_contigs.py"
