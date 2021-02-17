#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
from math import ceil

from haystac.workflow.scripts.utilities import PE, print_error

SUBSAMPLE_FIXED_READS = 200000


rule bowtie_align_db_single_end:
    input:
        fastq=(
            config["sample_output_dir"] + "/fastq_inputs/{read_mode}/{sample}_adRm.fastq.gz"
            if config["trim_adapters"]
            else config["fastq"]
        ),
        bt2idx=config["db_output"] + "/bowtie/chunk{chunk_num}.1.bt2l",
    log:
        config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_chunk{chunk_num}.log",
    output:
        bam_file=temp(config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_sorted_chunk{chunk_num}.bam"),
    params:
        index=config["db_output"] + "/bowtie/chunk{chunk_num}",
    threads: config["cores"]
    resources:
        mem_mb=(
            lambda wildcards: os.stat(config["db_output"] + f"/bowtie/chunk{wildcards.chunk_num}.fasta.gz").st_size * 5
        ),
    wildcard_constraints:
        read_mode="(SE|COLLAPSED)",
    message:
        "The filtering alignment for sample {wildcards.sample}, for index chunk number {wildcards.chunk_num} "
        "is being executed."
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -U {input.fastq} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"


rule bowtie_align_db_paired_end:
    input:
        fastq_r1=(
            config["sample_output_dir"] + "/fastq_inputs/PE/{sample}_R1_adRm.fastq.gz"
            if config["trim_adapters"]
            else config["fastq_r1"]
        ),
        fastq_r2=(
            config["sample_output_dir"] + "/fastq_inputs/PE/{sample}_R2_adRm.fastq.gz"
            if config["trim_adapters"]
            else config["fastq_r2"]
        ),
        bt2idx=config["db_output"] + "/bowtie/chunk{chunk_num}.1.bt2l",
    log:
        config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_chunk{chunk_num}.log",
    output:
        bam_file=temp(config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_sorted_chunk{chunk_num}.bam"),
    params:
        index=config["db_output"] + "/bowtie/chunk{chunk_num}",
    threads: config["cores"]
    resources:
        mem_mb=(
            lambda wildcards: os.stat(config["db_output"] + f"/bowtie/chunk{wildcards.chunk_num}.fasta.gz").st_size * 5
        ),
    wildcard_constraints:
        read_mode="(PE)",
    message:
        "The filtering alignment for sample {wildcards.sample}, for index chunk number {wildcards.chunk_num} "
        "is being executed."
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -1 {input.fastq_r1} -2 {input.fastq_r2} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"


def get_db_chunk_alignments(wildcards):
    """Get the sorted bam paths for each index chunk"""

    idx_chunk_total = ceil(float(open(config["db_output"] + "/bowtie/bt2_idx_chunk_num.txt").read().strip()))

    return expand(
        config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_sorted_chunk{chunk_num}.bam",
        read_mode=config["read_mode"],
        sample=wildcards.sample,
        chunk_num=range(1, idx_chunk_total + 1),
    )


rule merge_db_alignments:
    input:
        aln_path=get_db_chunk_alignments,
    log:
        config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_merge_bams.log",
    output:
        config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_sorted.bam",
    message:
        "Merging the bam files produced by the filtering alignment stage for sample {wildcards.sample}."
    wildcard_constraints:
        read_mode="[^-_]+",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools merge -f {output} {input.aln_path} 2> {log}"


rule extract_db_fastq_single_end:
    input:
        config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_sorted.bam",
    log:
        config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq.log",
    output:
        config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq.fastq.gz",
    message:
        "Extracting all the aligned reads for sample {wildcards.sample} and storing them in {output}. "
        "The log file can be found in {log}."
    conda:
        "../envs/extract_fastq.yaml"
    shell:
        "( samtools view -h -F 4 {input} | samtools fastq -c 6 - | seqkit rmdup -n -o {output} ) 2> {log}"


rule extract_db_fastq_paired_end:
    input:
        config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_sorted.bam",
    log:
        config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq.log",
    output:
        fastq_r1=config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq_R1.fastq.gz",
        fastq_r2=config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq_R2.fastq.gz",
        temp_r1=temp(config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_temp_R1_mapq.fastq.gz"),
        temp_r2=temp(config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_temp_R2_mapq.fastq.gz"),
    message:
        "Extracting the aligned reads for sample {wildcards.sample}."
    conda:
        "../envs/extract_fastq.yaml"
    shell:
        "( samtools view -h -F 4 {input} | "
        "  samtools fastq -c 6 -1 {output.temp_r1} -2 {output.temp_r2} -0 /dev/null -s /dev/null - &&"
        "  seqkit rmdup -n {output.temp_r1} -o {output.fastq_r1} &&"
        "  seqkit rmdup -n {output.temp_r2} -o {output.fastq_r2} "
        ") 2> {log}"


def get_extracted_db_fastq(wildcards):
    """Input function to assert that the fastqs are not empty"""

    if config["read_mode"] == PE:
        input_file = config["analysis_output_dir"] + f"/fastq/{wildcards.read_mode}/{wildcards.sample}_mapq_R1.fastq.gz"
    else:
        input_file = config["analysis_output_dir"] + f"/fastq/{wildcards.read_mode}/{wildcards.sample}_mapq.fastq.gz"

    if os.path.exists(input_file):
        if os.stat(input_file).st_size == 0:
            print_error(f"None of the reads in the sample file match any of the taxa in the database.")

    return input_file


rule average_fastq_read_len_single_end:
    input:
        get_extracted_db_fastq,
    log:
        config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq_readlen.log",
    output:
        config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq.readlen",
    message:
        "Calculating the average read length for sample {wildcards.sample}."
    conda:
        "../envs/seqtk.yaml"
    shell:
        "( seqtk sample {input} {SUBSAMPLE_FIXED_READS} | seqtk seq -A | grep -v '^>'"
        "| awk '{{count++; bases += length}} END {{print bases/count}}' 1> {output} ) 2> {log}"


rule average_fastq_read_len_paired_end:
    input:
        mate1=get_extracted_db_fastq,
        mate2=config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq_R2.fastq.gz",
    log:
        config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq_readlen.log",
    output:
        mate1=temp(config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq_R1.readlen"),
        mate2=temp(config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq_R2.readlen"),
        pair=config["analysis_output_dir"] + "/fastq/{read_mode}/{sample}_mapq_pair.readlen",
    message:
        "Calculating the average read length for sample {wildcards.sample}."
    conda:
        "../envs/seqtk.yaml"
    shell:
        "( seqtk sample {input.mate1} {SUBSAMPLE_FIXED_READS} | seqtk seq -A | grep -v '^>' "
        "| awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output.mate1} ) 2>> {log} && "
        "( seqtk sample {input.mate2} {SUBSAMPLE_FIXED_READS} | seqtk seq -A | grep -v '^>' "
        "| awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output.mate2} ) 2>> {log} && "
        "( cat {output.mate1} {output.mate2} "
        "| awk '{{sum += $1; n++ }} END {{if (n>0) print sum/n;}}' 1> {output.pair} ) 2>> {log} "
