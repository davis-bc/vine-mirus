rule taxonomy_checkv:
    """Extract viral contigs, assign NCBI taxonomy via DIAMOND, assess completeness with CheckV."""
    input:
        scores="results/intermediate/evidence/viral_scores.tsv",
        contigs=[get_contigs_path(s) for s in CONTIG_SETS],
        db=[] if config.get("test_mode") else [
            "databases/diamond/nr_viral.dmnd",
            "databases/taxdump/.done",
            "databases/checkv_db/.done",
        ],
    output:
        taxonomy="results/intermediate/taxonomy/viral_taxonomy.tsv",
        checkv="results/intermediate/checkv/quality_summary.tsv",
        viral_fa="results/viral_contigs.fa",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        diamond_db=config["diamond_refseq_viral_db"],
        taxdump=config["taxdump_dir"],
        checkv_db=config["checkv_db"],
        evalue=config["diamond_evalue"],
        max_target=config["diamond_max_target_seqs"],
        outfmt=config["diamond_outfmt"],
        checkv_outdir="results/intermediate/checkv",
        min_len=config["min_contig_length"],
        test_mode=config.get("test_mode", False),
    conda:
        "../../envs/viral_analysis.yaml"
    log:
        "logs/taxonomy_checkv.log"
    script:
        "../scripts/run_taxonomy_checkv.py"
