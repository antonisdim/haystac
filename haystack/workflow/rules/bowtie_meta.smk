#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

MIN_FRAG_LEN = 0
MAX_FRAG_LEN = 1000
META_ALN_MIN_SCORE_CONSTANT = -6

from haystack.workflow.scripts.utilities import get_total_paths, reads


def get_min_score(wildcards, input):
    """Get the min score dor the edit distance of the alignment."""
    return (
        round(float(open(input.readlen).read()) * float(config["mismatch_probability"])) * META_ALN_MIN_SCORE_CONSTANT
    )


rule align_taxon_single_end:
    input:
        fastq=config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.fastq.gz",
        db_idx=config["cache"] + "/ncbi/{orgname}/{accession}.1.bt2l",
        readlen=config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.readlen",
    log:
        config["analysis_output_dir"] + "/alignments/{sample}/SE/{orgname}/{accession}.log",
    output:
        bam_file=config["analysis_output_dir"] + "/alignments/{sample}/SE/{orgname}/{orgname}_{accession}.bam",
        bai_file=config["analysis_output_dir"] + "/alignments/{sample}/SE/{orgname}/{orgname}_{accession}.bam.bai",
    benchmark:
        repeat(
            "benchmarks/align_taxon_single_end_{sample}_{orgname}_{accession}.benchmark.txt", 1,
        )
    params:
        min_score=get_min_score,
        basename=config["cache"] + "/ncbi/{orgname}/{accession}",
    threads: config["bowtie2_threads"]
    message:
        "Aligning the filtered reads from sample {wildcards.sample} against taxon {wildcards.orgname}."
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "( bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
        "   --score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {threads} "
        "   -x {params.basename} -I {MIN_FRAG_LEN} -X {MAX_FRAG_LEN} -U {input.fastq} "
        "| samtools sort -O bam -o {output.bam_file} && samtools index {output.bam_file} "
        ") 2> {log}"


rule align_taxon_paired_end:
    input:
        fastq_r1=config["analysis_output_dir"] + "/fastq/PE/{sample}_R1_mapq.fastq.gz",
        fastq_r2=config["analysis_output_dir"] + "/fastq/PE/{sample}_R2_mapq.fastq.gz",
        db_idx=config["cache"] + "/ncbi/{orgname}/{accession}.1.bt2l",
        readlen=config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq_pair.readlen",
    log:
        config["analysis_output_dir"] + "/alignments/{sample}/PE/{orgname}/{accession}.log",
    output:
        bam_file=config["analysis_output_dir"] + "/alignments/{sample}/PE/{orgname}/{orgname}_{accession}.bam",
        bai_file=config["analysis_output_dir"] + "/alignments/{sample}/PE/{orgname}/{orgname}_{accession}.bam.bai",
    params:
        min_score=get_min_score,
        basename=config["cache"] + "/ncbi/{orgname}/{accession}",
    threads: config["bowtie2_threads"]
    message:
        "Aligning the filtered reads from sample {wildcards.sample} against taxon {wildcards.orgname}."
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "( bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
        "   --score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {threads} "
        "   -x {params.basename} -I {MIN_FRAG_LEN} -X {MAX_FRAG_LEN} -1 {input.fastq_r1} -2 {input.fastq_r2} "
        "| samtools sort -O bam -o {output.bam_file} && samtools index {output.bam_file} "
        ") 2> {log}"


rule all_alignments:
    input:
        [
            config["analysis_output_dir"]
            + "/alignments/{sample}/"
            + reads(config)
            + f"/{orgname}/{orgname}_{accession}.bam"
            for orgname, accession in get_total_paths(checkpoints, config)
        ],
    output:
        config["analysis_output_dir"] + "/alignments/{sample}_alignments.done",
    benchmark:
        repeat("benchmarks/all_alignments_{sample}.benchmark.txt", 1)
    message:
        "All metagenomic alignments are done."
    shell:
        "touch {output}"
