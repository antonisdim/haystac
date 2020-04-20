#!/usr/bin/env python
# -*- coding: utf-8 -*-


PAIRED_END_MODERN = True

SRA_LOOKUP = True
PE_ANCIENT = False
PE_MODERN = True
SE = False


##### Target rules #####


rule get_sra_fastq_se:
    output:
        temp("sra_data/{sample_accession}.fastq")
    log:
        temp("sra_data/{sample_accession}.log")
    params:
        extra=""
    threads: 6
    wrapper:
        "0.51.2/bio/sra-tools/fasterq-dump"



rule get_sra_fastq_pe:
    output:
        # the wildcard name must be accession, pointing to an SRA number
          temp("sra_data/{sample_accession}_1.fastq"),
          temp("sra_data/{sample_accession}_2.fastq")
    log:
        temp("sra_data/{sample_accession}.log")
    params:
        # optional extra arguments
        extra=""
    threads: 6  # defaults to 6
    wrapper:
        "0.51.2/bio/sra-tools/fasterq-dump"



rule compress_sra_fastq_se:
    input:
        "sra_data/{sample_accession}.fastq"
    log:
       "sra_data/{sample_accession}_compress.log"
    output:
        "sra_data/{sample_accession}.fastq.gz"
    shell:
        "gzip -c {input} > {output}"



rule compress_sra_fastq_pe:
    input:
        r1="sra_data/{sample_accession}_1.fastq",
        r2="sra_data/{sample_accession}_2.fastq"
    log:
        "sra_data/{sample_accession}_compress.log"
    output:
        r1="sra_data/{sample_accession}_R1.fastq.gz",
        r2="sra_data/{sample_accession}_R2.fastq.gz"
    shell:
        "gzip -c {input.r1} > {output.r1}; gzip -c {input.r2} > {output.r2}"



ruleorder: get_sra_fastq_pe > get_sra_fastq_se
ruleorder: compress_sra_fastq_pe > compress_sra_fastq_se


def get_inputs_for_adapterremoval_r1(wildcards):
    print(wildcards.sample_accession)

    if SRA_LOOKUP:
        if PE_ANCIENT or PE_MODERN:
            return "sra_data/{sample_accession}_R1.fastq.gz".format(sample_accession=wildcards.sample_accession)
        elif SE:
            return "sra_data/{sample_accession}.fastq.gz".format(sample_accession=wildcards.sample_accession)

    else:
        if PE_ANCIENT or PE_MODERN:
            return config['samples'][wildcards.sample_accession]['R1']
        elif SE:
            return config['samples'][wildcards.sample_accession]


def get_inputs_for_adapterremoval_r2(wildcards):

    if SRA_LOOKUP:
        if PE_ANCIENT or PE_MODERN:
            return "sra_data/{sample_accession}_R2.fastq.gz".format(sample_accession=wildcards.sample_accession)

    else:
        if PE_ANCIENT or PE_MODERN:
            return config['samples'][wildcards.sample_accession]['R2']



rule adapterremoval_single_end:
    input:
        fastq=get_inputs_for_adapterremoval_r1
    log:
        "fastq_inputs/{sample_accession}_adRm.log"
    output:
        "fastq_inputs/{sample_accession}_se_adRm.fastq.gz"
    benchmark:
        repeat("benchmarks/adapterremoval_single_end_{sample_accession}.benchmark.txt", 3)
    shell:
        "AdapterRemoval --file1 {input} --basename fastq_inputs/{wildcards.sample_accession} --gzip --minlength 15 "
        "--trimns; cat fastq_inputs/{wildcards.sample_accession}*truncated* > {output}"



rule adapterremoval_paired_end_ancient:
    input:
        fastq_r1=get_inputs_for_adapterremoval_r1,
        fastq_r2=get_inputs_for_adapterremoval_r2
    log:
        "fastq_inputs/{sample_accession}_adRm.log"
    output:
        "fastq_inputs/{sample_accession}_pe_adRm.fastq.gz"
    benchmark:
        repeat("benchmarks/adapterremoval_paired_end_ancient_{sample_accession}.benchmark.txt", 3)
    shell:
        "AdapterRemoval --file1 {input.fastq_r1}  --file2 {input.fastq_r2} "
        "--basename fastq_inputs/{wildcards.sample_accession} --gzip --collapse-deterministic  --minlength 15 "
        "--trimns; cat fastq_inputs/{wildcards.sample_accession}*collapsed* > {output}"



rule adapterremoval_paired_end_modern:
    input:
        fastq_r1=get_inputs_for_adapterremoval_r1,
        fastq_r2=get_inputs_for_adapterremoval_r2
    log:
        "fastq_inputs/{sample_accession}_adRm.log"
    output:
        fastq_r1="fastq_inputs/{sample_accession}_R1_adRm.fastq.gz",
        fastq_r2="fastq_inputs/{sample_accession}_R2_adRm.fastq.gz"
    benchmark:
        repeat("benchmarks/adapterremoval_paired_end_modern_{sample_accession}.benchmark.txt", 3)
    shell:
        "AdapterRemoval --file1 {input.fastq_r1}  --file2 {input.fastq_r2} "
        "--basename fastq_inputs/{wildcards.sample_accession} --gzip --minlength 15 --trimns; "
        "cat fastq_inputs/{wildcards.sample_accession}*pair1* > {output.fastq_r1}; "
        "cat fastq_inputs/{wildcards.sample_accession}*pair2* > {output.fastq_r2}"



# ruleorder: adapterremoval_paired_end_ancient > adapterremoval_single_end
# ruleorder: adapterremoval_paired_end_modern > adapterremoval_single_end

# rule choose_inputs:

    # todo input function in the bowtie.smk file with ifs for the returns - include the pair1 alignment options,
    #  constant static for any type of input, check paths or whether it is the config and then decide.
    #  Do that for adapter removal and bowtie.smk alignments
