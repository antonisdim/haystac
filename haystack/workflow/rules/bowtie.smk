#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from math import ceil

from haystack.workflow.scripts.utilities import (
    PE_MODERN,
    PE_ANCIENT,
    SE,
)

SUBSAMPLE_FIXED_READS = 200000


def get_inputs_for_bowtie_r1(wildcards):
    if config["trim_adapters"]:
        if config["read_mode"] == PE_MODERN:
            return config["sample_output_dir"] + f"/fastq_inputs/PE_mod/{wildcards.sample}_R1_adRm.fastq.gz"
        elif config["read_mode"] == PE_ANCIENT:
            return config["sample_output_dir"] + f"/fastq_inputs/PE_anc/{wildcards.sample}_adRm.fastq.gz"
        elif config["read_mode"] == SE:
            return config["sample_output_dir"] + f"/fastq_inputs/SE/{wildcards.sample}_adRm.fastq.gz"

    if config["read_mode"] == PE_MODERN:
        return config["fastq_r1"]
    elif config["read_mode"] == PE_ANCIENT or config["read_mode"] == SE:
        return config["fastq"]


def get_inputs_for_bowtie_r2(wildcards):
    if config["trim_adapters"]:
        if config["read_mode"] == PE_MODERN:
            return config["sample_output_dir"] + f"/fastq_inputs/PE_mod/{wildcards.sample}_R2_adRm.fastq.gz"

    if config["read_mode"] == PE_MODERN:
        return config["fastq_r2"]


rule bowtie_alignment_single_end:
    input:
        fastq=get_inputs_for_bowtie_r1,
        bt2idx=config["db_output"] + "/bowtie/chunk{chunk_num}.1.bt2l",
    log:
        config["analysis_output_dir"] + "/bam/{sample}_chunk{chunk_num}.log",
    output:
        bam_file=temp(config["analysis_output_dir"] + "/bam/SE_{sample}_sorted_chunk{chunk_num}.bam"),
    benchmark:
        repeat("benchmarks/bowtie_alignment_{sample}_chunk{chunk_num}.benchmark.txt", 1)
    params:
        index=config["db_output"] + "/bowtie/chunk{chunk_num}",
    threads: config["bowtie2_threads"]
    message:
        "The filtering alignment for sample {wildcards.sample}, for index chunk number {wildcards.chunk_num} "
        "is being executed."
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -U {input.fastq} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"


rule bowtie_alignment_paired_end:
    input:
        fastq_r1=get_inputs_for_bowtie_r1,
        fastq_r2=get_inputs_for_bowtie_r2,
        bt2idx=config["db_output"] + "/bowtie/chunk{chunk_num}.1.bt2l",
    log:
        config["analysis_output_dir"] + "/bam/{sample}_chunk{chunk_num}.log",
    output:
        bam_file=temp(config["analysis_output_dir"] + "/bam/PE_{sample}_sorted_chunk{chunk_num}.bam"),
    benchmark:
        repeat("benchmarks/bowtie_alignment_{sample}_chunk{chunk_num}.benchmark.txt", 1)
    params:
        index=config["db_output"] + "/bowtie/chunk{chunk_num}",
    threads: config["bowtie2_threads"]
    message:
        "The filtering alignment for sample {wildcards.sample}, for index chunk number {wildcards.chunk_num} "
        "is being executed."
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -1 {input.fastq_r1} -2 {input.fastq_r2} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"


def get_sorted_bam_paths(wildcards):
    """Get the sorted bam paths for each index chunk"""

    idx_chunk_total = ceil(float(open(config["db_output"] + "/bowtie/bt2_idx_chunk_num.txt").read().strip()))

    reads = ["PE"] if config["read_mode"] == PE_MODERN else ["SE"]

    return expand(
        config["analysis_output_dir"] + "/bam/{reads}_{sample}_sorted_chunk{chunk_num}.bam",
        reads=reads,
        sample=wildcards.sample,
        chunk_num=range(1, idx_chunk_total + 1),
    )


rule merge_bams:
    input:
        aln_path=get_sorted_bam_paths,
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


rule extract_fastq_single_end:
    input:
        config["analysis_output_dir"] + "/bam/SE_{sample}_sorted.bam",
    log:
        config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.log",
    output:
        config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.fastq.gz",
    benchmark:
        repeat("benchmarks/extract_fastq_single_end_{sample}.benchmark.txt", 1)
    message:
        "Extracting all the aligned reads for sample {wildcards.sample} and storing them in {output}. "
        "The log file can be found in {log}."
    conda:
        "../envs/extract_fastq.yaml"
    shell:
        "( samtools view -h -F 4 {input} | samtools fastq -c 6 - | seqkit rmdup -n -o {output} ) 2> {log}"


rule extract_fastq_paired_end:
    input:
        config["analysis_output_dir"] + "/bam/PE_{sample}_sorted.bam",
    log:
        config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq.log",
    output:
        fastq_r1=config["analysis_output_dir"] + "/fastq/PE/{sample}_R1_mapq.fastq.gz",
        fastq_r2=config["analysis_output_dir"] + "/fastq/PE/{sample}_R2_mapq.fastq.gz",
        temp_r1=temp(config["analysis_output_dir"] + "/fastq/PE/{sample}_temp_R1_mapq.fastq.gz"),
        temp_r2=temp(config["analysis_output_dir"] + "/fastq/PE/{sample}_temp_R2_mapq.fastq.gz"),
    benchmark:
        repeat("benchmarks/extract_fastq_paired_end_{sample}.benchmark.txt", 1)
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


rule average_fastq_read_len_single_end:
    input:
        config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.fastq.gz",
    log:
        config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq_readlen.log",
    output:
        config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.readlen",
    benchmark:
        repeat("benchmarks/average_fastq_read_len_single_end_{sample}.benchmark.txt", 1)
    message:
        "Calculating the average read length for sample {wildcards.sample}."
    conda:
        "../envs/seqtk.yaml"
    shell:
        "( seqtk sample {input} {SUBSAMPLE_FIXED_READS} | seqtk seq -A | grep -v '^>'"
        "| awk '{{count++; bases += length}} END {{print bases/count}}' 1> {output} ) 2> {log}"


rule average_fastq_read_len_paired_end:
    input:
        mate1=config["analysis_output_dir"] + "/fastq/PE/{sample}_R1_mapq.fastq.gz",
        mate2=config["analysis_output_dir"] + "/fastq/PE/{sample}_R2_mapq.fastq.gz",
    log:
        config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq_readlen.log",
    output:
        mate1=temp(config["analysis_output_dir"] + "/fastq/{sample}_R1_mapq.readlen"),
        mate2=temp(config["analysis_output_dir"] + "/fastq/{sample}_R2_mapq.readlen"),
        pair=config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq_pair.readlen",
    benchmark:
        repeat("benchmarks/average_fastq_read_len_paired_end_{sample}.benchmark.txt", 1)
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
