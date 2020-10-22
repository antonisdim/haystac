#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


def get_inputs_for_count_fastq_len(wildcards):
    if config["trim_adapters"]:
        if config["PE_MODERN"]:
            return config["sample_output_dir"] + f"/fastq_inputs/PE_mod/{wildcards.sample}_R1_adRm.fastq.gz"
        elif config["PE_ANCIENT"]:
            return config["sample_output_dir"] + f"/fastq_inputs/PE_anc/{wildcards.sample}_adRm.fastq.gz"
        elif config["SE"]:
            return config["sample_output_dir"] + f"/fastq_inputs/SE/{wildcards.sample}_adRm.fastq.gz"

    else:
        if config["PE_MODERN"]:
            return config["fastq_R1"]
        elif config["PE_ANCIENT"]:
            return config["fastq"]
        elif config["SE"]:
            return config["fastq"]


rule count_fastq_length:
    input:
        fastq=get_inputs_for_count_fastq_len,
    log:
        config["sample_output_dir"] + "/fastq_inputs/meta/{sample}.log",
    output:
        config["sample_output_dir"] + "/fastq_inputs/meta/{sample}.size",
    benchmark:
        repeat("benchmarks/count_fastq_length_{sample}.benchmark.txt", 1)
    message:
        "Counting the number of reads in sample {wildcards.sample}."
    conda:
        "../envs/seqtk.yaml"
    shell:
        "(seqtk seq -A {input.fastq} | grep -v '^>' | wc -l 1> {output} ) 2> {log}"
