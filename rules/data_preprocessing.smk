#!/usr/bin/env python
# -*- coding: utf-8 -*-


SRA_LOOKUP = config["SRA_LOOKUP"]
PE_ANCIENT = config["PE_ANCIENT"]
PE_MODERN = config["PE_MODERN"]
SE = config["SE"]

##### Target rules #####

# TODO move all SRA rules into a new `sra.smk` file, and rename this file `adapterremoval.smk`
rule get_sra_fastq_se:
    output:
        temp("sra_data/SE/{accession}.fastq"),
    log:
        temp("sra_data/SE/{accession}.log"),
    params:
        extra="",
    threads: 6
    message:
        "Download SRA file {output} for accession {wildcards.accession}. The log file can be found in {log}."
    wrapper:
        "0.51.2/bio/sra-tools/fasterq-dump"  # TODO we're not using wrappers anywhere else... be consistent!


rule get_sra_fastq_pe:
    output:
        temp("sra_data/PE/{accession}_1.fastq"), # the wildcard name must be accession, pointing to an SRA number
        temp("sra_data/PE/{accession}_2.fastq"),
    log:
        temp("sra_data/PE/{accession}.log"),
    params:
        extra="", # optional extra arguments
    threads: 6 # defaults to 6
    message:
        "Download SRA files {output} for accession {wildcards.accession}. The log file can be found in {log}."
    wrapper:
        "0.51.2/bio/sra-tools/fasterq-dump"  # TODO we're not using wrappers anywhere else... be consistent!


rule compress_sra_fastq_se:
    input:
        "sra_data/SE/{accession}.fastq",
    log:
        "sra_data/SE/{accession}_compress.log",  # TODO log file is not informative
    output:
        "sra_data/SE/{accession}.fastq.gz",
    message:
        "Compressing the raw fastq file {input} for accession {wildcards.accession} and storing it in {output}. "
        "The log file can be found in {log}."
    shell:
        "gzip -c {input} 1> {output} 2> {log}"


rule compress_sra_fastq_pe:
    input:
        r1="sra_data/PE/{accession}_1.fastq",
        r2="sra_data/PE/{accession}_2.fastq",
    log:
        "sra_data/PE/{accession}_compress.log",  # TODO not only is this log file not informative but you immediately overwrite it then you compress r2
    output:
        r1="sra_data/PE/{accession}_R1.fastq.gz",
        r2="sra_data/PE/{accession}_R2.fastq.gz",
    message:
        "Compressing the raw fastq files {input.r1} and {input.r2} for accession {wildcards.accession} "
        "and storing them in {output.r1} and {output.r2}. "
        "The log file can be found in {log}."
    shell:
        "gzip -c {input.r1} 1> {output.r1} 2> {log}; "
        "gzip -c {input.r2} 1> {output.r2} 2> {log}"


# TODO delete commented out code
# ruleorder: get_sra_fastq_pe > get_sra_fastq_se
# ruleorder: compress_sra_fastq_pe > compress_sra_fastq_se


def get_inputs_for_adapterremoval_r1(wildcards):
    # print(wildcards.accession)

    # TODO there is no need for this input function... you can set config["sample_fastq_R1"] = "sra_data/PE/{accession}_R1.fastq.gz" when in SRA mode
    if SRA_LOOKUP:
        if PE_ANCIENT or PE_MODERN:
            return "sra_data/PE/{accession}_R1.fastq.gz".format(
                accession=wildcards.accession
            )
        elif SE:
            return "sra_data/SE/{accession}.fastq.gz".format(
                accession=wildcards.accession
            )

    else:
        if PE_ANCIENT or PE_MODERN:
            return config["sample_fastq_R1"]
        elif SE:
            return config["sample_fastq"]


def get_inputs_for_adapterremoval_r2(wildcards):
    if SRA_LOOKUP:
        if PE_ANCIENT or PE_MODERN:
            return "sra_data/PE/{accession}_R2.fastq.gz".format(
                accession=wildcards.accession
            )
        else:
            return ""

    else:
        if PE_ANCIENT or PE_MODERN:
            return config["sample_fastq_R2"]
        else:
            return ""


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
        "--trimns; cat fastq_inputs/SE/{wildcards.accession}*truncated* 1> {output} 2> {log}"  # TODO NEVER USE FILE GLOBBING!!!!
        # TODO what if someone runs a sample with an accession that is one letter longer than an existing accession? e.g. 'SRR1031289' and 'SRR10312891'


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
        "--trimns; cat fastq_inputs/PE_anc/{wildcards.accession}*collapsed* 1> {output} 2> {log}" # TODO NEVER USE FILE GLOBBING!!!!


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
        "cat fastq_inputs/PE_mod/{wildcards.accession}*pair1* 1> {output.fastq_r1} 2> {log}; " # TODO NEVER USE FILE GLOBBING!!!!
        "cat fastq_inputs/PE_mod/{wildcards.accession}*pair2* 1> {output.fastq_r2} 2> {log}" # TODO NEVER USE FILE GLOBBING!!!!


# ruleorder: adapterremoval_paired_end_modern > adapterremoval_paired_end_ancient
