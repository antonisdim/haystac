#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from haystac.workflow.scripts.utilities import PE


def get_inputs_for_count_fastq_len(wildcards):
    if config["trim_adapters"]:
        if config["read_mode"] == PE:
            return (
                config["sample_output_dir"] + f"/fastq_inputs/{config['read_mode']}/{wildcards.sample}_R1_adRm.fastq.gz"
            )
        else:
            return config["sample_output_dir"] + f"/fastq_inputs/{config['read_mode']}/{wildcards.sample}_adRm.fastq.gz"

    else:
        return config["fastq_r1"] or config["fastq"]


rule count_fastq_length:
    input:
        fastq=get_inputs_for_count_fastq_len,
    output:
        config["sample_output_dir"] + "/fastq_inputs/meta/{sample}.size",
    message:
        "Counting the number of reads in sample {wildcards.sample}."
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk seq -A {input.fastq} | grep -v '^>' | wc -l 1> {output}"
