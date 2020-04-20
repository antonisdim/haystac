#!/usr/bin/env python
# -*- coding: utf-8 -*-


SRA_LOOKUP = True
PE_ANCIENT = True
PE_MODERN = False
SE = False


##### Target rules #####


rule get_sra_fastq_se:
    output:
        temp("sra_data/SE/{accession}.fastq")
    log:
        temp("sra_data/SE/{accession}.log")
    params:
        extra=""
    threads: 6
    wrapper:
        "0.51.2/bio/sra-tools/fasterq-dump"



rule get_sra_fastq_pe:
    output:
        # the wildcard name must be accession, pointing to an SRA number
          temp("sra_data/PE/{accession}_1.fastq"),
          temp("sra_data/PE/{accession}_2.fastq")
    log:
        temp("sra_data/PE/{accession}.log")
    params:
        # optional extra arguments
        extra=""
    threads: 6  # defaults to 6
    wrapper:
        "0.51.2/bio/sra-tools/fasterq-dump"



rule compress_sra_fastq_se:
    input:
        "sra_data/SE/{accession}.fastq"
    log:
       "sra_data/SE/{accession}_compress.log"
    output:
        "sra_data/SE/{accession}.fastq.gz"
    shell:
        "gzip -c {input} 1> {output} 2> {log}"



rule compress_sra_fastq_pe:
    input:
        r1="sra_data/PE/{accession}_1.fastq",
        r2="sra_data/PE/{accession}_2.fastq"
    log:
        "sra_data/PE/{accession}_compress.log"
    output:
        r1="sra_data/PE/{accession}_R1.fastq.gz",
        r2="sra_data/PE/{accession}_R2.fastq.gz"
    shell:
        "gzip -c {input.r1} 1> {output.r1} 2> {log}; gzip -c {input.r2} 1> {output.r2} 2> {log}"



ruleorder: get_sra_fastq_pe > get_sra_fastq_se
ruleorder: compress_sra_fastq_pe > compress_sra_fastq_se


def get_inputs_for_adapterremoval_r1(wildcards):
    print(wildcards.accession)

    if SRA_LOOKUP:
        if PE_ANCIENT or PE_MODERN:
            return "sra_data/PE/{accession}_R1.fastq.gz".format(accession=wildcards.accession)
        elif SE:
            return "sra_data/SE/{accession}.fastq.gz".format(accession=wildcards.accession)

    else:
        if PE_ANCIENT or PE_MODERN:
            return config['samples'][wildcards.accession]['R1']
        elif SE:
            return config['samples'][wildcards.accession]


def get_inputs_for_adapterremoval_r2(wildcards):

    if SRA_LOOKUP:
        if PE_ANCIENT or PE_MODERN:
            return "sra_data/PE/{accession}_R2.fastq.gz".format(accession=wildcards.accession)

    else:
        if PE_ANCIENT or PE_MODERN:
            return config['samples'][wildcards.accession]['R2']



rule adapterremoval_single_end:
    input:
        fastq=get_inputs_for_adapterremoval_r1
    log:
        "fastq_inputs/SE/{accession}_adRm.log"
    output:
        "fastq_inputs/SE/{accession}_adRm.fastq.gz"
    benchmark:
        repeat("benchmarks/adapterremoval_single_end_{accession}.benchmark.txt", 3)
    shell:
        "AdapterRemoval --file1 {input} --basename fastq_inputs/SE/{wildcards.accession} --gzip --minlength 15 "
        "--trimns; cat fastq_inputs/SE/{wildcards.accession}*truncated* 1> {output} 2> {log}"



rule adapterremoval_paired_end_ancient:
    input:
        fastq_r1=get_inputs_for_adapterremoval_r1,
        fastq_r2=get_inputs_for_adapterremoval_r2
    log:
        "fastq_inputs/PE/{accession}_adRm.log"
    output:
        "fastq_inputs/PE/{accession}_adRm.fastq.gz"
    benchmark:
        repeat("benchmarks/adapterremoval_paired_end_ancient_{accession}.benchmark.txt", 3)
    shell:
        "AdapterRemoval --file1 {input.fastq_r1}  --file2 {input.fastq_r2} "
        "--basename fastq_inputs/PE/{wildcards.accession} --gzip --collapse-deterministic  --minlength 15 "
        "--trimns; cat fastq_inputs/PE/{wildcards.accession}*collapsed* 1> {output} 2> {log}"



rule adapterremoval_paired_end_modern:
    input:
        fastq_r1=get_inputs_for_adapterremoval_r1,
        fastq_r2=get_inputs_for_adapterremoval_r2
    log:
        "fastq_inputs/PE/{accession}_adRm.log"
    output:
        fastq_r1="fastq_inputs/PE/{accession}_R1_adRm.fastq.gz",
        fastq_r2="fastq_inputs/PE/{accession}_R2_adRm.fastq.gz"
    benchmark:
        repeat("benchmarks/adapterremoval_paired_end_modern_{accession}.benchmark.txt", 3)
    shell:
        "AdapterRemoval --file1 {input.fastq_r1}  --file2 {input.fastq_r2} "
        "--basename fastq_inputs/PE/{wildcards.accession} --gzip --minlength 15 --trimns; "
        "cat fastq_inputs/PE/{wildcards.accession}*pair1* 1> {output.fastq_r1} 2> {log}; "
        "cat fastq_inputs/PE/{wildcards.accession}*pair2* 1> {output.fastq_r2} 2> {log}"



ruleorder: adapterremoval_paired_end_modern > adapterremoval_paired_end_ancient


# # rule choose_inputs:
#
#     # todo input function in the bowtie.smk file with ifs for the returns - include the pair1 alignment options,
#     #  constant static for any type of input, check paths or whether it is the config and then decide.
#     #  Do that for adapter removal and bowtie.smk alignments
