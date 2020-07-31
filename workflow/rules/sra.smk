#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

MESSAGE_SUFFIX = "(output: {output} and log: {log})" if config["debug"] else ""

##### Target rules #####


rule get_sra_fastq_se:
    output:
        temp(config["sample_output_dir"] + "/sra_data/SE/{accession}.fastq"),
    log:
        temp(config["sample_output_dir"] + "/sra_data/SE/{accession}.log"),
    params:
        extra="",
    threads: 6
    message:
        "Download SRA file for accession {wildcards.accession} {MESSAGE_SUFFIX}"
    conda:
        "../envs/sra_tools.yaml"
    shell:
        "fasterq-dump --split-files {wildcards.accession} --outdir " + config["sample_output_dir"] + "/sra_data/SE/",


rule get_sra_fastq_pe:
    output:
        temp(config["sample_output_dir"] + "/sra_data/PE/{accession}_1.fastq"), # the wildcard name must be accession, pointing to an SRA number
        temp(config["sample_output_dir"] + "/sra_data/PE/{accession}_2.fastq"),
    log:
        temp(config["sample_output_dir"] + "/sra_data/PE/{accession}.log"),
    params:
        extra="", # optional extra arguments
    threads: 6 # defaults to 6
    message:
        "Download SRA files for accession {wildcards.accession} {MESSAGE_SUFFIX}"
    conda:
        "../envs/sra_tools.yaml"
    shell:
        "fasterq-dump --split-files {wildcards.accession} --outdir " + config["sample_output_dir"] + "/sra_data/PE/"


rule compress_sra_fastq_se:
    input:
        config["sample_output_dir"] + "/sra_data/SE/{accession}.fastq",
    output:
        config["sample_output_dir"] + "/sra_data/SE/{accession}.fastq.gz",
    message:
        "Compressing the raw fastq file {input} {MESSAGE_SUFFIX}"
    shell:
        "gzip -c {input} 1> {output}"


rule compress_sra_fastq_pe:
    input:
        r1=config["sample_output_dir"] + "/sra_data/PE/{accession}_1.fastq",
        r2=config["sample_output_dir"] + "/sra_data/PE/{accession}_2.fastq",
    output:
        r1=config["sample_output_dir"] + "/sra_data/PE/{accession}_R1.fastq.gz",
        r2=config["sample_output_dir"] + "/sra_data/PE/{accession}_R2.fastq.gz",
    message:
        "Compressing the raw fastq files {input.r1} and {input.r2} {MESSAGE_SUFFIX}"
    shell:
        "gzip -c {input.r1} 1> {output.r1}; "
        "gzip -c {input.r2} 1> {output.r2}"
