rule evidence_score:
    input:
        genomad=expand("results/samples/{sample}/genomad/{sample}_contigs_summary/{sample}_contigs_virus_summary.tsv", sample=CONTIG_SETS),
        virsorter2=expand("results/samples/{sample}/virsorter2/final-viral-score.tsv", sample=CONTIG_SETS),
        pprmeta=expand("results/samples/{sample}/pprmeta/{sample}_pprmeta.csv", sample=CONTIG_SETS),
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
