#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


rule get_sra_fastq_se:
    output:
        temp(config["sample_output_dir"] + "/sra_data/SE/{accession}.fastq"),
    log:
        temp(config["sample_output_dir"] + "/sra_data/SE/{accession}.log"),
    threads: 6
    message:
        "Download SRA file for accession {wildcards.accession}."
    conda:
        "../envs/sra_tools.yaml"
    params:
        basename=config["sample_output_dir"] + "/sra_data/SE/",
    shell:
        # TODO this downloads a massive SRA file to `~/ncbi/public/sra/{accession}.sra.cache`
        "fasterq-dump {wildcards.accession}"
        " --split-files"
        " --threads {threads}"
        " --temp {params.basename}"
        " --outdir {params.basename} &> {log}"


rule get_sra_fastq_pe:
    output:
        temp(config["sample_output_dir"] + "/sra_data/PE/{accession}_1.fastq"),
        temp(config["sample_output_dir"] + "/sra_data/PE/{accession}_2.fastq"),
    log:
        temp(config["sample_output_dir"] + "/sra_data/PE/{accession}.log"),
    threads: 6
    message:
        "Download SRA files for accession {wildcards.accession}."
    conda:
        "../envs/sra_tools.yaml"
    params:
        basename=config["sample_output_dir"] + "/sra_data/PE/",
    shell:
        # TODO this downloads a massive SRA file to `~/ncbi/public/sra/{accession}.sra.cache`
        "fasterq-dump {wildcards.accession}"
        " --split-files"
        " --threads {threads}"
        " --temp {params.basename}"
        " --outdir {params.basename} &> {log}"


rule compress_sra_fastq_se:
    input:
        config["sample_output_dir"] + "/sra_data/SE/{accession}.fastq",
    output:
        config["sample_output_dir"] + "/sra_data/SE/{accession}.fastq.gz",
    message:
        "Compressing the raw fastq file {input}."
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
        "Compressing the raw fastq files {input.r1} and {input.r2}."
    shell:
        "gzip -c {input.r1} 1> {output.r1}; "
        "gzip -c {input.r2} 1> {output.r2}"
