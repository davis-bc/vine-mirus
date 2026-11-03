import os

# databases.smk — download and build every tool database under databases/.
#
# Each rule produces a sentinel (.done) or the real artifact, and tool rules
# declare those as inputs (empty in test mode). Snakemake's touch() writes the
# sentinel only after the shell block succeeds, so every rule VERIFIES its
# download before claiming success — several of these downloaders exit 0 having
# written nothing usable, or having written it one directory deeper than their
# own tools then look for it.
#
# The nesting each downloader actually produces (confirmed on Jordan's runs):
#   genomad download-database X  -> X/genomad_db/          (not X)
#   checkv download_database X   -> X/checkv-db-v1.N/      (version in the name)
#   download-db.sh X             -> X/databases/           (+ a sibling files/)
#   virsorter setup              -> the CURRENT DIRECTORY, ignoring --db-dir
# config/config.yaml therefore points some keys at the nested path, and the
# CheckV rule resolves its versioned directory at run time.
#
# Pre-fetch everything without running the pipeline:
#   snakemake --snakefile workflow/Snakefile databases --cores 4 --use-conda
#
# NOTE: these rules download tens to hundreds of GB from NCBI and the tool
# authors' servers. They cannot be exercised by the test-mode run, so treat the
# first real execution of each as the thing being tested, one at a time.

# ── Aggregate target ──────────────────────────────────────────────────────────
def _database_targets():
    """Sentinels for the databases the enabled tools actually need."""
    targets = []
    if TOOLS.get("genomad"):
        targets.append("databases/genomad_db/.done")
    if TOOLS.get("virsorter2"):
        targets.append("databases/virsorter2_db/.done")
    if TOOLS.get("vibrant"):
        targets.append("databases/vibrant_db/.done")
    if TOOLS.get("cenote_taker3"):
        targets.append("databases/cenote_taker3_db/.done")
    if TOOLS.get("diamond_refseq_viral"):
        targets.append(config["diamond_refseq_viral_db"])
    if TOOLS.get("diamond_nr"):
        targets += [config["diamond_nr_db"], "databases/taxdump/.done"]
    if TOOLS.get("pprmeta"):
        # Literal rather than PPRMETA_SCRIPT: pprmeta.smk is included after
        # this file, and this function runs at include time.
        targets.append("databases/pprmeta/PPR_Meta.py")
    # CheckV runs on the candidates regardless of which detection tools are on.
    targets.append("databases/checkv_db/.done")
    return targets


rule databases:
    """Download/build every database the enabled tools need, then stop."""
    input:
        _database_targets(),


# ── geNomad ───────────────────────────────────────────────────────────────────
rule download_genomad_db:
    """geNomad database (~3 GB). Creates a nested genomad_db/ subdirectory."""
    output:
        done=touch("databases/genomad_db/.done"),
    params:
        outdir=lambda wc, output: os.path.dirname(output.done),
        # The real database directory, and what config.genomad_db points at.
        db=config["genomad_db"],
    threads: 1
    resources:
        mem_mb=4000,
        runtime=240,
    conda:
        "../../envs/genomad.yaml"
    log:
        "logs/databases/genomad_db.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        genomad download-database {params.outdir} > {log} 2>&1
        # geNomad's own end-to-end run reads version.txt from the nested dir;
        # pointing a tool at the parent fails with FileNotFoundError.
        if [ ! -s "{params.db}/version.txt" ]; then
            echo "ERROR: expected {params.db}/version.txt after download." >> {log}
            echo "Inspect {params.outdir} - geNomad's layout may have changed." >> {log}
            exit 1
        fi
        echo "geNomad database verified at {params.db}" >> {log}
        """


# ── VirSorter2 ────────────────────────────────────────────────────────────────
rule download_virsorter2_db:
    """VirSorter2 database (~10 GB) plus a one-off build of its internal env.

    Two confirmed problems are worked around here:

    1. VirSorter2's setup.smk sets `db_dir = os.path.abspath(os.getcwd())` and
       never reads --db-dir / --config, so the database lands wherever the
       command was invoked from while still reporting success. The only
       mechanism it respects is the working directory, so this rule cd's into
       the target first.
    2. VirSorter2 launches its own nested Snakemake with --use-conda at RUN
       time. Letting that happen inside a pipeline job means every parallel
       chunk races to build the same env, which fails under Slurm. The smoke
       run below builds it once, here, serially. It is allowed to fail: a tiny
       synthetic contig legitimately yields no viral calls, and all we need is
       for the environment to have been created.
    """
    output:
        done=touch("databases/virsorter2_db/.done"),
    params:
        outdir=lambda wc, output: os.path.dirname(output.done),
        groups=",".join(config["virsorter2_groups"]),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=600,
    conda:
        "../../envs/virsorter2.yaml"
    log:
        "logs/databases/virsorter2_db.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        log=$(realpath {log})
        outdir=$(realpath {params.outdir})

        # cd first: setup.smk uses the cwd, not --db-dir (see docstring).
        ( cd "$outdir" && virsorter setup --db-dir . --jobs {threads} ) > "$log" 2>&1

        if [ -z "$(find "$outdir" -name '*.hmm' -print -quit)" ]; then
            echo "ERROR: no HMM files under $outdir after setup." >> "$log"
            echo "The download reported success but wrote nothing usable." >> "$log"
            exit 1
        fi

        # Pre-build VirSorter2's nested conda env with a throwaway run.
        smoke=$(mktemp -d)
        trap 'rm -rf "$smoke"' EXIT
        {{
            echo ">smoke_contig_1"
            for i in $(seq 1 60); do
                echo "ATGCGTACGTAGCTAGCTAGGCTAACGTTACGGATCCGGATTACAGGCATGCATCGATCGA"
            done
        }} > "$smoke/smoke.fa"
        echo "Pre-building VirSorter2's internal conda env via a smoke run" >> "$log"
        virsorter run -i "$smoke/smoke.fa" -w "$smoke/work" \
            --include-groups {params.groups} -j 1 --db-dir "$outdir" all \
            >> "$log" 2>&1 || echo "[NOTE] smoke run did not complete; see above" >> "$log"

        echo "VirSorter2 database verified at $outdir" >> "$log"
        """


# ── CheckV ────────────────────────────────────────────────────────────────────
rule download_checkv_db:
    """CheckV database (~3 GB). Creates a versioned checkv-db-v1.N/ inside."""
    output:
        done=touch("databases/checkv_db/.done"),
    params:
        outdir=lambda wc, output: os.path.dirname(output.done),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=240,
    conda:
        "../../envs/checkv.yaml"
    log:
        "logs/databases/checkv_db.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        checkv download_database {params.outdir} > {log} 2>&1
        # The version is in the directory name, so checkv.smk resolves it with
        # find at run time rather than pinning a version here.
        db=$(find {params.outdir} -maxdepth 1 -name 'checkv-db-*' -type d | head -1)
        if [ -z "$db" ] || [ ! -s "$db/genome_db/checkv_reps.dmnd" ]; then
            echo "ERROR: no checkv-db-* directory with a built genome_db under" >> {log}
            echo "{params.outdir}. Inspect it manually." >> {log}
            exit 1
        fi
        echo "CheckV database verified at $db" >> {log}
        """


# ── VIBRANT ───────────────────────────────────────────────────────────────────
rule download_vibrant_db:
    """VIBRANT database (~11 GB).

    VIBRANT_run.py needs TWO paths: -d <db>/databases (the HMM sets it
    downloads) and -m <db>/files (small annotation tables that ship inside the
    package). The downloader is `download-db.sh` — hyphen, singular "db"; the
    name in VIBRANT's own README does not exist — and it writes into a
    `databases` subdirectory of the path given, not the path itself. Passing
    the parent is what produces "could not identify KEGG HMM files".

    `files/` ships with the installed package rather than being downloaded, so
    if the downloader did not place it here it is copied out of the conda env.
    Normalizing the layout in this rule keeps vibrant.smk's paths simple.
    """
    output:
        done=touch("databases/vibrant_db/.done"),
    params:
        outdir=lambda wc, output: os.path.dirname(output.done),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=600,
    conda:
        "../../envs/vibrant.yaml"
    log:
        "logs/databases/vibrant_db.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        outdir=$(realpath {params.outdir})
        log=$(realpath {log})

        download-db.sh "$outdir" > "$log" 2>&1

        # Some builds download in place inside the package instead; adopt
        # whichever "databases" directory actually has the HMM sets.
        if [ ! -s "$outdir/databases/KEGG_profiles_prokaryotes.HMM" ]; then
            src=$(find "$CONDA_PREFIX" -type d -name databases -path '*VIBRANT*' \
                  -exec test -s '{{}}/KEGG_profiles_prokaryotes.HMM' \; -print -quit || true)
            if [ -n "$src" ]; then
                echo "[NOTE] adopting VIBRANT databases from $src" >> "$log"
                mkdir -p "$outdir/databases"
                cp -rn "$src/." "$outdir/databases/"
            fi
        fi

        if [ ! -d "$outdir/files" ]; then
            src=$(find "$CONDA_PREFIX" -type d -name files -path '*VIBRANT*' -print -quit || true)
            if [ -n "$src" ]; then
                echo "[NOTE] copying VIBRANT files/ from $src" >> "$log"
                cp -r "$src" "$outdir/files"
            fi
        fi

        for required in "$outdir/databases/KEGG_profiles_prokaryotes.HMM" "$outdir/files"; do
            if [ ! -e "$required" ]; then
                echo "ERROR: missing $required after setup." >> "$log"
                echo "VIBRANT needs -d <db>/databases and -m <db>/files; one of" >> "$log"
                echo "them was not produced. Inspect $outdir." >> "$log"
                exit 1
            fi
        done
        echo "VIBRANT database verified at $outdir" >> "$log"
        """


# ── Cenote-Taker3 ─────────────────────────────────────────────────────────────
rule download_cenote_taker3_db:
    """Cenote-Taker3 standard databases (~3 GB decompressed).

    The downloader is `get_ct3_dbs`, an installed entry point. The three
    hhsuite databases (--hhCDD 6 GB, --hhPFAM 5 GB, --hhPDB 56 GB) are
    optional and skipped: they would dominate the disk budget and are not used
    by the hallmark-domain detection this pipeline relies on.

    Cenote-Taker3 also reads a CENOTE_DBS environment variable. Rather than
    `conda env config vars set` — which only takes effect on the next
    activation, and silently reverts when the env is rebuilt — cenote_taker3.smk
    exports it in the job itself alongside --cenote-dbs.
    """
    output:
        done=touch("databases/cenote_taker3_db/.done"),
    params:
        outdir=lambda wc, output: os.path.dirname(output.done),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=480,
    conda:
        "../../envs/cenote_taker3.yaml"
    log:
        "logs/databases/cenote_taker3_db.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        get_ct3_dbs -o {params.outdir} \
            --hmm T --hallmark_tax T --refseq_tax T --mmseqs_cdd T --domain_list T \
            > {log} 2>&1
        if [ -z "$(find {params.outdir} -type f -print -quit)" ]; then
            echo "ERROR: get_ct3_dbs left {params.outdir} empty." >> {log}
            exit 1
        fi
        du -sh {params.outdir} >> {log} 2>&1 || true
        echo "Cenote-Taker3 databases verified at {params.outdir}" >> {log}
        """


# ── NCBI taxonomy ─────────────────────────────────────────────────────────────
rule download_taxdump:
    """NCBI taxonomy dump: nodes.dmp + names.dmp (~70 MB compressed).

    Used twice: taxopy resolves DIAMOND hits to lineages with it, and
    `diamond makedb --taxonnodes/--taxonnames` needs the same two files.
    """
    output:
        done=touch("databases/taxdump/.done"),
        nodes="databases/taxdump/nodes.dmp",
        names="databases/taxdump/names.dmp",
    params:
        outdir=lambda wc, output: os.path.dirname(output.nodes),
        url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
    conda:
        "../../envs/download.yaml"
    log:
        "logs/databases/taxdump.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        curl -fsSL {params.url} -o {params.outdir}/taxdump.tar.gz 2> {log}
        tar -xzf {params.outdir}/taxdump.tar.gz -C {params.outdir} nodes.dmp names.dmp 2>> {log}
        rm -f {params.outdir}/taxdump.tar.gz
        test -s {output.nodes} && test -s {output.names}
        echo "taxdump verified: $(wc -l < {output.nodes}) nodes" >> {log}
        """


rule download_prot_accession2taxid:
    """protein accession -> taxid map, for `diamond makedb --taxonmap` (~2 GB).

    Supplying this at makedb time is what puts real staxids on every hit.
    Jordan's blastx-based route joined accessions to taxids with an awk filter
    afterwards and produced staxids=0 for every hit — a well-formed result in
    which nothing could ever be classified as viral.
    """
    output:
        taxonmap="databases/taxonomy/prot.accession2taxid.gz",
    params:
        url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=240,
    conda:
        "../../envs/download.yaml"
    log:
        "logs/databases/prot_accession2taxid.log"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.taxonmap})
        curl -fsSL {params.url} -o {output.taxonmap} 2> {log}
        # Truncated transfers are the common failure; gzip -t catches them.
        gzip -t {output.taxonmap} 2>> {log}
        echo "prot.accession2taxid verified ($(du -h {output.taxonmap} | cut -f1))" >> {log}
        """


# ── DIAMOND: RefSeq viral proteins ────────────────────────────────────────────
rule build_diamond_refseq_viral_db:
    """DIAMOND database of NCBI RefSeq viral proteins, taxonomy-enabled (~2 GB).

    RefSeq viral is split across viral.1, viral.2, ... .protein.faa.gz files,
    and the number of parts grows between releases. Jordan's setup fetched only
    viral.1, which silently truncated the reference: every search then ran
    against a partial database and simply found less. This loop takes parts
    until one 404s, and fails if it got none.
    """
    input:
        taxonmap="databases/taxonomy/prot.accession2taxid.gz",
        nodes="databases/taxdump/nodes.dmp",
        names="databases/taxdump/names.dmp",
    output:
        dmnd=config["diamond_refseq_viral_db"],
    params:
        stem=lambda wc, output: os.path.splitext(output.dmnd)[0],
        workdir=lambda wc, output: os.path.join(os.path.dirname(output.dmnd), "refseq_viral_parts"),
        ftp="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=480,
    conda:
        "../../envs/diamond.yaml"
    log:
        "logs/databases/diamond_refseq_viral.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.workdir}
        n=0
        for i in $(seq 1 200); do
            f="viral.${{i}}.protein.faa.gz"
            if curl -fsS -o "{params.workdir}/$f" "{params.ftp}/$f" 2>> {log}; then
                echo "downloaded $f" >> {log}
                n=$((n+1))
            else
                rm -f "{params.workdir}/$f"
                break
            fi
        done
        if [ "$n" -eq 0 ]; then
            echo "ERROR: no RefSeq viral protein parts downloaded from {params.ftp}" >> {log}
            exit 1
        fi
        echo "RefSeq viral protein parts: $n" >> {log}

        zcat {params.workdir}/viral.*.protein.faa.gz | \
            diamond makedb --in - --db {params.stem} \
                --taxonmap {input.taxonmap} \
                --taxonnodes {input.nodes} --taxonnames {input.names} \
                --threads {threads} >> {log} 2>&1
        rm -rf {params.workdir}
        test -s {output.dmnd}
        echo "DIAMOND RefSeq viral database built at {output.dmnd}" >> {log}
        """


# ── DIAMOND: nr ───────────────────────────────────────────────────────────────
# nr is ~250 GB compressed and takes many hours to download and build, so it is
# never fetched implicitly: set download_nr: true to opt in, or point
# diamond_nr_db at a copy that already exists (e.g. a shared cluster path).
if DOWNLOAD_NR:

    rule build_diamond_nr_db:
        """DIAMOND database of NCBI nr, taxonomy-enabled. Hours, hundreds of GB."""
        input:
            taxonmap="databases/taxonomy/prot.accession2taxid.gz",
            nodes="databases/taxdump/nodes.dmp",
            names="databases/taxdump/names.dmp",
        output:
            dmnd=config["diamond_nr_db"],
        params:
            stem=lambda wc, output: os.path.splitext(output.dmnd)[0],
            fasta=lambda wc, output: os.path.join(os.path.dirname(output.dmnd), "nr.gz"),
            url="https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz",
        threads: 16
        resources:
            mem_mb=64000,
            runtime=2880,
        conda:
            "../../envs/diamond.yaml"
        log:
            "logs/databases/diamond_nr.log"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.dmnd})
            # -C - resumes a partial transfer: re-fetching nr from scratch
            # after a makedb failure would cost hours for nothing.
            if [ ! -s {params.fasta} ]; then
                curl -fsSL -C - {params.url} -o {params.fasta} 2> {log}
            else
                echo "nr.gz already present; not re-downloading" > {log}
            fi
            gzip -t {params.fasta} 2>> {log}
            diamond makedb --in {params.fasta} --db {params.stem} \
                --taxonmap {input.taxonmap} \
                --taxonnodes {input.nodes} --taxonnames {input.names} \
                --threads {threads} >> {log} 2>&1
            test -s {output.dmnd}
            echo "DIAMOND nr database built at {output.dmnd}" >> {log}
            """

else:

    rule diamond_nr_db_missing:
        """Explain how to provide nr instead of silently starting a 250 GB download."""
        output:
            dmnd=config["diamond_nr_db"],
        conda:
            "../../envs/download.yaml"
        log:
            "logs/databases/diamond_nr.log"
        shell:
            """
            cat >&2 <<MSG
{output.dmnd} does not exist and download_nr is false.

The nr search is Tier B's cellular-false-positive check. To provide it:
  - set download_nr: true in config/config.yaml to build it here
    (hundreds of GB and many hours), or
  - point diamond_nr_db at an existing taxonomy-enabled nr.dmnd, or
  - set virus_identification.tools.diamond_nr: false to skip the search.
MSG
            exit 1
            """
