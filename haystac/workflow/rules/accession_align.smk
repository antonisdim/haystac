#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
from math import ceil

MIN_FRAG_LEN = 0
MAX_FRAG_LEN = 1000
META_ALN_MIN_SCORE_CONSTANT = -6  # bowtie2 default base mismatch penalty value

from haystac.workflow.scripts.utilities import get_final_db_paths


# noinspection PyUnusedLocal,PyShadowingBuiltins,PyShadowingNames
def get_min_score(wildcards, input):
    """Get the min score dor the edit distance of the alignment."""
    return ceil(float(open(input.readlen).read()) * float(config["mismatch_probability"])) * META_ALN_MIN_SCORE_CONSTANT


rule bowtie_align_accession_single_end:
    input:
        fastq=config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq.fastq.gz",
        db_idx=config["cache"] + "/ncbi/{orgname}/{accession}.1.bt2l",
        readlen=config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq.readlen",
    log:
        config["analysis_output_dir"] + "/alignments/{sample}/{read_mode}/{orgname}/{accession}.log",
    output:
        bam_file=(
            config["analysis_output_dir"] + "/alignments/{sample}/{read_mode}/{orgname}/{orgname}_{accession}.bam"
        ),
        bai_file=(
            config["analysis_output_dir"] + "/alignments/{sample}/{read_mode}/{orgname}/{orgname}_{accession}.bam.bai"
        ),
    params:
        min_score=get_min_score,
        basename=config["cache"] + "/ncbi/{orgname}/{accession}",
    threads: config["bowtie2_threads_aln"]
    resources:
        mem_mb=(
            lambda wildcards: os.stat(
                config["cache"] + f"/ncbi/{wildcards.orgname}/{wildcards.accession}.fasta.gz"
            ).st_size
            * 5
        ),
    wildcard_constraints:
        read_mode="(SE|COLLAPSED)",
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


rule bowtie_align_accession_paired_end:
    input:
        fastq_r1=config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq_R1.fastq.gz",
        fastq_r2=config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq_R2.fastq.gz",
        db_idx=config["cache"] + "/ncbi/{orgname}/{accession}.1.bt2l",
        readlen=config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq_pair.readlen",
    log:
        config["analysis_output_dir"] + "/alignments/{sample}/{read_mode}/{orgname}/{accession}.log",
    output:
        bam_file=(
            config["analysis_output_dir"] + "/alignments/{sample}/{read_mode}/{orgname}/{orgname}_{accession}.bam"
        ),
        bai_file=(
            config["analysis_output_dir"] + "/alignments/{sample}/{read_mode}/{orgname}/{orgname}_{accession}.bam.bai"
        ),
    params:
        min_score=get_min_score,
        basename=config["cache"] + "/ncbi/{orgname}/{accession}",
    threads: config["bowtie2_threads_aln"]
    resources:
        mem_mb=(
            lambda wildcards: os.stat(
                config["cache"] + f"/ncbi/{wildcards.orgname}/{wildcards.accession}.fasta.gz"
            ).st_size
            * 5
        ),
    wildcard_constraints:
        read_mode="(PE)",
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


def get_accession_alignments(_):
    """Get paths for individual bam files"""
    return [
        config["analysis_output_dir"]
        + "/alignments/{sample}/"
        + config["read_mode"]
        + f"/{orgname}/{orgname}_{accession}.bam"
        for orgname, accession in get_final_db_paths(checkpoints)
    ]


rule align_all_accessions:
    input:
        get_accession_alignments,
    output:
        config["analysis_output_dir"] + "/alignments/{sample}_{read_mode}_alignments.done",
    message:
        "All metagenomic alignments are done."
    shell:
        "touch {output}"
