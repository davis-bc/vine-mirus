import os

# preprocess.smk — read QC for every sample with reads.
#
# qc_reads runs fastp → [host removal] → [dedup] inside ONE job: these always
# run back-to-back on the same reads, so separate rules would only add
# scheduling overhead and full-size intermediate FASTQs on disk.
# Output reads (clean_reads()) feed both assembly and abundance mapping;
# normalization for assembly happens in assembly.smk and never touches them.

_PRE = config["preprocessing"]
HOST_REMOVAL = bool(_PRE["host_removal"]["enabled"])
HOST_REFS = list(_PRE["host_removal"].get("reference_genomes") or [])
DEDUP = bool(_PRE["deduplication"]["enabled"])

if HOST_REMOVAL and not HOST_REFS:
    raise ValueError("preprocessing.host_removal.enabled is true but reference_genomes is empty")

HOST_INDEX = "results/host/host_sr.mmi"


if HOST_REMOVAL:
    rule host_index:
        """Concatenate host reference(s) and build a minimap2 short-read index.

        minimap2 is used instead of bwa: indexing a ~1.5 Gb fish genome takes
        minutes instead of ~1 h, and mapping is faster at similar sensitivity
        for read-level host filtering.
        """
        input:
            refs=HOST_REFS,
        output:
            index=HOST_INDEX,
        threads: 8
        resources:
            mem_mb=32000,
            runtime=120,
        conda:
            "../../envs/qc_assembly.yaml"
        log:
            "logs/host_index.log"
        shell:
            """
            tmp=$(mktemp -d results/host/tmp.XXXXXX)
            cat {input.refs} > $tmp/host.fa
            minimap2 -x sr -t {threads} -d {output.index} $tmp/host.fa 2> {log}
            rm -rf $tmp
            """


rule qc_reads:
    input:
        r1=get_r1,
        r2=get_r2,
        host_index=[HOST_INDEX] if HOST_REMOVAL else [],
    output:
        r1="results/samples/{sample}/qc/{sample}_R1.fastq.gz",
        r2="results/samples/{sample}/qc/{sample}_R2.fastq.gz",
        json="results/samples/{sample}/qc/{sample}_fastp.json",
        html="results/samples/{sample}/qc/{sample}_fastp.html",
    threads: 8
    resources:
        # minimap2 holds the whole host index in memory; clumpify is RAM-hungry
        mem_mb=32000 if (HOST_REMOVAL or DEDUP) else 8000,
        runtime=120 if (HOST_REMOVAL or DEDUP) else 60,
    params:
        min_len=_PRE["fastp"]["min_length"],
        qual=_PRE["fastp"]["qualified_quality"],
        host=int(HOST_REMOVAL),
        dedup=int(DEDUP),
        tmpdir="results/samples/{sample}/qc/tmp",
    conda:
        "../../envs/qc_assembly.yaml"
    log:
        "logs/qc_reads/{sample}.log"
    shell:
        """
        set -euo pipefail
        rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
        cur1={params.tmpdir}/fastp_R1.fastq.gz
        cur2={params.tmpdir}/fastp_R2.fastq.gz

        fastp \
            -i {input.r1} -I {input.r2} \
            -o $cur1 -O $cur2 \
            --json {output.json} --html {output.html} \
            --length_required {params.min_len} \
            --qualified_quality_phred {params.qual} \
            --thread {threads} \
            2> {log}

        if [ {params.host} -eq 1 ]; then
            echo "[host removal] minimap2 vs {input.host_index}" >> {log}
            minimap2 -ax sr -t {threads} {input.host_index} $cur1 $cur2 2>> {log} \
                | samtools fastq -f 12 -F 256 -n -c 1 \
                    -1 {params.tmpdir}/nohost_R1.fastq.gz -2 {params.tmpdir}/nohost_R2.fastq.gz \
                    -0 /dev/null -s /dev/null - 2>> {log}
            cur1={params.tmpdir}/nohost_R1.fastq.gz
            cur2={params.tmpdir}/nohost_R2.fastq.gz
        fi

        if [ {params.dedup} -eq 1 ]; then
            echo "[dedup] clumpify" >> {log}
            # Java pre-commits its heap: cap at 80% of min(requested, available)
            avail=$(awk '/MemAvailable/ {{print int($2/1024)}}' /proc/meminfo)
            req={resources.mem_mb}; xmx=$(( (req < avail ? req : avail) * 8 / 10 ))
            clumpify.sh in1=$cur1 in2=$cur2 \
                out1={params.tmpdir}/dedup_R1.fastq.gz out2={params.tmpdir}/dedup_R2.fastq.gz \
                dedupe=t threads={threads} -Xmx${{xmx}}m 2>> {log}
            cur1={params.tmpdir}/dedup_R1.fastq.gz
            cur2={params.tmpdir}/dedup_R2.fastq.gz
        fi

        mv $cur1 {output.r1}
        mv $cur2 {output.r2}
        rm -rf {params.tmpdir}

        # Exit 0 is not proof of success: refuse to emit empty read files.
        # (zcat || true): head closing the pipe SIGPIPEs zcat (exit 141) under pipefail
        n=$( (zcat {output.r1} || true) | head -c 1 | wc -c)
        if [ "$n" -eq 0 ]; then
            echo "ERROR: no reads left for {wildcards.sample} after QC/host/dedup" >> {log}
            exit 1
        fi
        """


rule multiqc:
    input:
        fastp=expand("results/samples/{sample}/qc/{sample}_fastp.json", sample=READS_SAMPLES),
        quast=expand("results/samples/{unit}/assembly/quast/report.tsv", unit=ASSEMBLY_UNITS)
            if config["assembly"].get("quast", False) else [],
    output:
        html="results/qc/multiqc_report.html",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        outdir=lambda wc, output: os.path.dirname(output.html),
    conda:
        "../../envs/reports.yaml"
    log:
        "logs/multiqc.log"
    shell:
        "multiqc --force -o {params.outdir} -n multiqc_report.html {input} 2> {log}"
