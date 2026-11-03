# assembly.smk — one assembly per "assembly unit" (see Snakefile):
#   coassembly mode : unit "coasm" = all ASSEMBLE_SAMPLES pooled
#   per_sample mode : unit = each sample
#
# Normalized reads are temp() and used for assembly ONLY. Abundance mapping
# (coverm.smk) always uses the un-normalized clean_reads().

import re

_ASM = config["assembly"]
_NORM = config["preprocessing"]["normalization"]
NORMALIZE = bool(_NORM["enabled"])
ASSEMBLER = _ASM["assembler"]


# Only real assembly units may be produced by these rules; otherwise a
# user-provided contigs path under results/ could be "re-assembled".
_UNIT_RE = "|".join(re.escape(u) for u in ASSEMBLY_UNITS) or "^$"


def _unit_clean_reads(mate):
    def f(wildcards):
        return [clean_reads(s, mate) for s in unit_samples(wildcards.unit)]
    return f


if NORMALIZE:
    rule normalize_reads:
        """bbnorm the reads of one assembly unit (pooled first for coasm)."""
        input:
            r1=_unit_clean_reads(1),
            r2=_unit_clean_reads(2),
        wildcard_constraints:
            unit=_UNIT_RE,
        output:
            r1=temp("results/samples/{unit}/assembly/norm/norm_R1.fastq.gz"),
            r2=temp("results/samples/{unit}/assembly/norm/norm_R2.fastq.gz"),
        threads: 8
        resources:
            mem_mb=32000,
            runtime=120,
        params:
            target=_NORM["target_depth"],
            min_depth=_NORM["min_depth"],
            tmpdir="results/samples/{unit}/assembly/norm/tmp",
        conda:
            "../../envs/qc_assembly.yaml"
        log:
            "logs/normalize_reads/{unit}.log"
        shell:
            """
            set -euo pipefail
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
            # bbnorm reads its input twice (count pass + normalize pass), so a
            # pool of several samples is concatenated to a temp file first.
            # gzip members concatenate into a valid gzip stream.
            set -- {input.r1}
            if [ $# -gt 1 ]; then
                cat {input.r1} > {params.tmpdir}/pool_R1.fastq.gz
                cat {input.r2} > {params.tmpdir}/pool_R2.fastq.gz
                in1={params.tmpdir}/pool_R1.fastq.gz; in2={params.tmpdir}/pool_R2.fastq.gz
            else
                in1={input.r1}; in2={input.r2}
            fi
            # Java pre-commits its heap: cap at 80% of min(requested, available)
            avail=$(awk '/MemAvailable/ {{print int($2/1024)}}' /proc/meminfo)
            req={resources.mem_mb}; xmx=$(( (req < avail ? req : avail) * 8 / 10 ))
            bbnorm.sh in1=$in1 in2=$in2 out1={output.r1} out2={output.r2} \
                target={params.target} min={params.min_depth} \
                threads={threads} -Xmx${{xmx}}m 2> {log}
            rm -rf {params.tmpdir}
            """


def _assembly_reads(mate):
    if NORMALIZE:
        return lambda wc: f"results/samples/{wc.unit}/assembly/norm/norm_R{mate}.fastq.gz"
    return _unit_clean_reads(mate)


rule assemble:
    wildcard_constraints:
        unit=_UNIT_RE,
    input:
        r1=_assembly_reads(1),
        r2=_assembly_reads(2),
    output:
        contigs="results/samples/{unit}/assembly/{unit}_contigs.fa",
    threads: 16
    resources:
        mem_mb=64000,
        runtime=240,
    params:
        assembler=ASSEMBLER,
        min_len=_ASM["min_contig_len"],
        outdir="results/samples/{unit}/assembly/tmp",
        r1_csv=lambda wc, input: ",".join(input.r1 if isinstance(input.r1, list) else [input.r1]),
        r2_csv=lambda wc, input: ",".join(input.r2 if isinstance(input.r2, list) else [input.r2]),
        mem_gb=lambda wc, resources: max(1, resources.mem_mb // 1000),
    conda:
        "../../envs/qc_assembly.yaml"
    log:
        "logs/assemble/{unit}.log"
    shell:
        """
        set -euo pipefail
        rm -rf {params.outdir}
        if [ {params.assembler} = megahit ]; then
            # MEGAHIT takes comma-separated libraries, so no pooling copy is needed
            megahit -1 {params.r1_csv} -2 {params.r2_csv} \
                --min-contig-len {params.min_len} \
                -t {threads} -m {resources.mem_mb}000000 \
                -o {params.outdir} 2> {log}
            mv {params.outdir}/final.contigs.fa {output.contigs}
        else
            # metaSPAdes accepts only one paired library
            mkdir -p {params.outdir}
            set -- {input.r1}
            if [ $# -gt 1 ]; then
                cat {input.r1} > {params.outdir}/pool_R1.fastq.gz
                cat {input.r2} > {params.outdir}/pool_R2.fastq.gz
                in1={params.outdir}/pool_R1.fastq.gz; in2={params.outdir}/pool_R2.fastq.gz
            else
                in1={input.r1}; in2={input.r2}
            fi
            metaspades.py -1 $in1 -2 $in2 -t {threads} -m {params.mem_gb} \
                -o {params.outdir}/spades 2> {log}
            seqkit seq -m {params.min_len} {params.outdir}/spades/contigs.fasta \
                > {output.contigs} 2>> {log}
        fi
        rm -rf {params.outdir}
        """


if _ASM.get("quast", False):
    rule quast:
        input:
            contigs="results/samples/{unit}/assembly/{unit}_contigs.fa",
        output:
            report="results/samples/{unit}/assembly/quast/report.tsv",
        threads: 4
        resources:
            mem_mb=8000,
            runtime=60,
        params:
            outdir="results/samples/{unit}/assembly/quast",
            min_len=_ASM["min_contig_len"],
        conda:
            "../../envs/reports.yaml"
        log:
            "logs/quast/{unit}.log"
        shell:
            "quast.py {input.contigs} -o {params.outdir} -t {threads} "
            "--min-contig {params.min_len} --labels {wildcards.unit} 2> {log}"
