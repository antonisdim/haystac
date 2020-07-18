#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


##### Target rules #####


def get_inputs_for_adapterremoval_r1(wildcards):
    # print(wildcards.accession)

    # TODO there is no need for this input function... you can set config["sample_fastq_R1"] = "sra_data/PE/{accession}_R1.fastq.gz" when in SRA mode
    if config["SRA_LOOKUP"]:
        if config["PE_ANCIENT"] or config["PE_MODERN"]:
            return "sra_data/PE/{accession}_R1.fastq.gz".format(
                accession=wildcards.accession
            )
        elif config["SE"]:
            return "sra_data/SE/{accession}.fastq.gz".format(
                accession=wildcards.accession
            )

    else:
        if config["PE_ANCIENT"] or config["PE_MODERN"]:
            return config["sample_fastq_R1"]
        elif config["SE"]:
            return config["sample_fastq"]


def get_inputs_for_adapterremoval_r2(wildcards):
    if config["SRA_LOOKUP"]:
        if config["PE_ANCIENT"] or config["PE_MODERN"]:
            return "sra_data/PE/{accession}_R2.fastq.gz".format(
                accession=wildcards.accession
            )

    else:
        if config["PE_ANCIENT"] or config["PE_MODERN"]:
            return config["sample_fastq_R2"]


rule adapterremoval_single_end:
    input:
        fastq=get_inputs_for_adapterremoval_r1,
    log:
        "fastq_inputs/SE/{accession}_adRm.log",
    output:
        "fastq_inputs/SE/{accession}_adRm.fastq.gz",
    benchmark:
        repeat("benchmarks/adapterremoval_single_end_{accession}.benchmark.txt", 1)
    message:
        "Trimming sequencing adapters from file {input.fastq}, using AdapterRemoval. "
        "The trimmed reads can be found in {output}, and the "
        "log file can be found in {log}."
    shell:
        "AdapterRemoval --file1 {input} --basename fastq_inputs/SE/{wildcards.accession} --gzip --minlength 15 "
        "--trimns; cat fastq_inputs/SE/{wildcards.accession}.truncated.gz 1> {output} 2> {log}"


rule adapterremoval_paired_end_ancient:
    input:
        fastq_r1=get_inputs_for_adapterremoval_r1,
        fastq_r2=get_inputs_for_adapterremoval_r2,
    log:
        "fastq_inputs/PE_anc/{accession}_adRm.log",
    output:
        "fastq_inputs/PE_anc/{accession}_adRm.fastq.gz",
    benchmark:
        repeat("benchmarks/adapterremoval_paired_end_ancient_{accession}.benchmark.txt", 1)
    message:
        "Trimming sequencing adapters and collapsing reads from files {input.fastq_r1} and {input.fastq_r2}, "
        "using AdapterRemoval. "
        "The trimmed reads can be found in {output}, and the "
        "log file can be found in {log}."
    shell:
        "AdapterRemoval --file1 {input.fastq_r1}  --file2 {input.fastq_r2} "
        "--basename fastq_inputs/PE_anc/{wildcards.accession} --gzip --collapse-deterministic  --minlength 15 "
        "--trimns; cat fastq_inputs/PE_anc/{wildcards.accession}.collapsed.gz "
        "fastq_inputs/PE_anc/{wildcards.accession}.collapsed.truncated.gz 1> {output} 2> {log}"


rule adapterremoval_paired_end_modern:
    input:
        fastq_r1=get_inputs_for_adapterremoval_r1,
        fastq_r2=get_inputs_for_adapterremoval_r2,
    log:
        "fastq_inputs/PE_mod/{accession}_adRm.log",
    output:
        fastq_r1="fastq_inputs/PE_mod/{accession}_R1_adRm.fastq.gz",
        fastq_r2="fastq_inputs/PE_mod/{accession}_R2_adRm.fastq.gz",
    benchmark:
        repeat("benchmarks/adapterremoval_paired_end_modern_{accession}.benchmark.txt", 1)
    message:
        "Trimming sequencing adapters from files {input.fastq_r1} and {input.fastq_r2}, "
        "using AdapterRemoval. "
        "The trimmed reads can be found in {output.fastq_r1} and {output.fastq_r2}, and the "
        "log file can be found in {log}."
    shell:
        "AdapterRemoval --file1 {input.fastq_r1}  --file2 {input.fastq_r2} "
        "--basename fastq_inputs/PE_mod/{wildcards.accession} --gzip --minlength 15 --trimns; "
        "cat fastq_inputs/PE_mod/{wildcards.accession}.pair1.truncated.gz 1> {output.fastq_r1} 2> {log}; "
        "cat fastq_inputs/PE_mod/{wildcards.accession}.pair2.truncated.gz 1> {output.fastq_r2} 2> {log}"
