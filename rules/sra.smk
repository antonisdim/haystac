#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


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
        "Download SRA file {output} for accession {wildcards.accession}. The log file can be found in {log}."
    wrapper:
        "0.51.2/bio/sra-tools/fasterq-dump" # TODO we're not using wrappers anywhere else... be consistent!


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
        "Download SRA files {output} for accession {wildcards.accession}. The log file can be found in {log}."
    wrapper:
        "0.51.2/bio/sra-tools/fasterq-dump" # TODO we're not using wrappers anywhere else... be consistent!


rule compress_sra_fastq_se:
    input:
        config["sample_output_dir"] + "/sra_data/SE/{accession}.fastq",
    output:
        config["sample_output_dir"] + "/sra_data/SE/{accession}.fastq.gz",
    message:
        "Compressing the raw fastq file {input} for accession {wildcards.accession} and storing it in {output}."
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
        "Compressing the raw fastq files {input.r1} and {input.r2} for accession {wildcards.accession} "
        "and storing them in {output.r1} and {output.r2}."
    shell:
        "gzip -c {input.r1} 1> {output.r1}; "
        "gzip -c {input.r2} 1> {output.r2}"
