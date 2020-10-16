#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

MESSAGE_SUFFIX = "(output: {output} and log: {log})" if config["debug"] else ""


rule adapterremoval_single_end:
    input:
        fastq=config["fastq"] if config["fastq"] else config["fastq_r1"],
    log:
        config["sample_output_dir"] + "/fastq_inputs/SE/{accession}_adRm.log",
    output:
        config["sample_output_dir"] + "/fastq_inputs/SE/{accession}_adRm.fastq.gz",
    benchmark:
        repeat("benchmarks/adapterremoval_single_end_{accession}.benchmark.txt", 1)
    message:
        "Trimming sequencing adapters from file {input.fastq} {MESSAGE_SUFFIX}"
    conda:
        "../envs/adapterremoval.yaml"
    shell:
        "(AdapterRemoval --file1 {input} --basename {config[sample_output_dir]}/fastq_inputs/SE/{wildcards.accession} "
        "--gzip --minlength 15 --trimns; "
        "cat {config[sample_output_dir]}/fastq_inputs/SE/{wildcards.accession}.truncated.gz 1> {output}) 2> {log}"


rule adapterremoval_paired_end_ancient:
    input:
        fastq_r1=config["fastq_r1"] if config["fastq_r1"] else config["fastq"],
        fastq_r2=config["fastq_r2"],
    log:
        config["sample_output_dir"] + "/fastq_inputs/PE_anc/{accession}_adRm.log",
    output:
        config["sample_output_dir"] + "/fastq_inputs/PE_anc/{accession}_adRm.fastq.gz",
    benchmark:
        repeat("benchmarks/adapterremoval_paired_end_ancient_{accession}.benchmark.txt", 1)
    message:
        "Trimming sequencing adapters and collapsing reads from files {input.fastq_r1} and {input.fastq_r2} "
        "{MESSAGE_SUFFIX}"
    conda:
        "../envs/adapterremoval.yaml"
    shell:
        "(AdapterRemoval --file1 {input.fastq_r1}  --file2 {input.fastq_r2} "
        "--basename {config[sample_output_dir]}/fastq_inputs/PE_anc/{wildcards.accession} --gzip "
        "--collapse-deterministic  --minlength 15 "
        "--trimns; cat {config[sample_output_dir]}/fastq_inputs/PE_anc/{wildcards.accession}.collapsed.gz "+
        "{config[sample_output_dir]}/fastq_inputs/PE_anc/{wildcards.accession}.collapsed.truncated.gz 1> {output}) "
        "2> {log}"


rule adapterremoval_paired_end_modern:
    input:
        fastq_r1=config["fastq_r1"] if config["fastq_r1"] else config["fastq"],
        fastq_r2=config["fastq_r2"],
    log:
        config["sample_output_dir"] + "/fastq_inputs/PE_mod/{accession}_adRm.log",
    output:
        fastq_r1=config["sample_output_dir"] + "/fastq_inputs/PE_mod/{accession}_R1_adRm.fastq.gz",
        fastq_r2=config["sample_output_dir"] + "/fastq_inputs/PE_mod/{accession}_R2_adRm.fastq.gz",
    benchmark:
        repeat("benchmarks/adapterremoval_paired_end_modern_{accession}.benchmark.txt", 1)
    message:
        "Trimming sequencing adapters from files {input.fastq_r1} and {input.fastq_r2} "
        "{MESSAGE_SUFFIX}"
    conda:
        "../envs/adapterremoval.yaml"
    shell:
        "(AdapterRemoval --file1 {input.fastq_r1}  --file2 {input.fastq_r2} "
        "--basename {config[sample_output_dir]}/fastq_inputs/PE_mod/{wildcards.accession} --gzip --minlength 15 --trimns; "
        "cat {config[sample_output_dir]}/fastq_inputs/PE_mod/{wildcards.accession}.pair1.truncated.gz "
        "1> {output.fastq_r1}; "
        "cat {config[sample_output_dir]}/fastq_inputs/PE_mod/{wildcards.accession}.pair2.truncated.gz "
        "1> {output.fastq_r2}) 2> {log}"
