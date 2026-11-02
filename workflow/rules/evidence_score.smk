rule evidence_score:
    input:
        genomad=tool_tables("genomad"),
        virsorter2=tool_tables("virsorter2"),
        pprmeta=tool_tables("pprmeta"),
        # Consumed by the evidence rewrite; listed here so the tools run and
        # their outputs are produced/validated.
        vibrant=tool_tables("vibrant"),
        cenote_taker3=tool_tables("cenote_taker3"),
        diamond_refseq_viral=tool_tables("diamond_refseq_viral"),
    output:
        scores="results/intermediate/evidence/viral_scores.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        samples=CONTIG_SETS,
        threshold=config["viral_score_threshold"],
        weights=config["tool_weights"],
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/evidence_score.log"
    script:
        "../scripts/score_evidence.py"
