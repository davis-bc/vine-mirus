import os

# Tier A collation: every enabled detection tool -> one per-contig table, plus
# the candidate FASTA that the expensive Tier B steps run on.

rule evidence_score:
    input:
        # Declared so Snakemake reruns when the shared parsing library changes
        # (it is imported at runtime via params.scripts, not tracked otherwise).
        lib=f"{SCRIPTS}/evidence_lib.py",
        contigs=[get_contigs_path(s) for s in CONTIG_SETS],
        genomad=tool_tables("genomad"),
        virsorter2=tool_tables("virsorter2"),
        vibrant=tool_tables("vibrant"),
        cenote_taker3=tool_tables("cenote_taker3"),
        pprmeta=tool_tables("pprmeta"),
        diamond_refseq_viral=tool_tables("diamond_refseq_viral"),
    output:
        scores="results/intermediate/evidence/evidence_tierA.tsv",
        candidates="results/intermediate/evidence/candidates.fa",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        # Derived from the tracked lib input so there is one source of truth
        scripts=lambda wc, input: os.path.dirname(input.lib),
        tools=TOOLS,
        thresholds={
            "genomad_min_score": config["virus_identification"]["genomad_min_score"],
            "virsorter2_min_score": config["virus_identification"]["virsorter2_min_score"],
            "pprmeta_min_score": config["virus_identification"]["pprmeta_min_score"],
        },
        diamond=config["virus_identification"]["diamond"],
        weights=config["tool_weights"],
        collation=config["collation"],
        include_phage=INCLUDE_PHAGE,
        strict_checks=STRICT_CHECKS,
        threshold=config["viral_score_threshold"],
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/evidence_score.log"
    script:
        "../scripts/score_evidence.py"
