import os

# Tier B collation: add the nr search and CheckV, assign taxonomy, finalize the
# host-domain and phage calls, and write the reviewer-facing tables.

rule evidence_final:
    input:
        # Declared so Snakemake reruns when the shared parsing library changes
        # (it is imported at runtime via params.scripts, not tracked otherwise).
        lib=f"{SCRIPTS}/evidence_lib.py",
        scores="results/intermediate/evidence/evidence_tierA.tsv",
        candidates="results/intermediate/evidence/candidates.fa",
        checkv="results/intermediate/checkv/quality_summary.tsv",
        diamond_nr=tool_tables("diamond_nr"),
        taxdump=[] if TEST_MODE or not TOOLS.get("diamond_nr") else ["databases/taxdump/.done"],
    output:
        evidence="results/evidence_table.tsv",
        phage="results/phage_evidence_table.tsv",
        viral_fa="results/viral_contigs.fa",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        # Derived from the tracked lib input so there is one source of truth
        scripts=lambda wc, input: os.path.dirname(input.lib),
        test_mode=TEST_MODE,
        taxdump="" if TEST_MODE else config["taxdump_dir"],
        diamond=config["virus_identification"]["diamond"],
        confidence_tiers=config["collation"]["confidence_tiers"],
        strict_checks=STRICT_CHECKS,
        include_phage=INCLUDE_PHAGE,
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/evidence_final.log"
    script:
        "../scripts/evidence_final.py"
