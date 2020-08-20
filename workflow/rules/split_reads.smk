#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from scripts.rip_utilities import get_total_paths, normalise_name

MESSAGE_SUFFIX = "(output: {output} and log: {log})" if config["debug"] else ""


##### Target rules #####


rule get_dirichlet_reads:
    input:
        bam_file=config["analysis_output_dir"] + "/alignments/{sample}/{reads}/{orgname}/{orgname}_{accession}.bam",
        dirichlet_matrix=config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
    output:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/{orgname}/{orgname}_{accession}_dirichlet_{reads}.bam",
    log:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/{orgname}/{orgname}_{accession}_dirichlet_{reads}.log",
    benchmark:
        repeat(
            "benchmarks/get_dirichlet_reads_{sample}_{reads}_{orgname}_{accession}.benchmark.txt", 1,
        )
    message:
        "Preparing bam files with the Dirichlet assigned reads for taxon {wildcards.orgname} "
        "for sample {wildcards.sample} {MESSAGE_SUFFIX}"
    conda:
        "../envs/dirichlet_reads.yaml"
    script:
        "../scripts/get_dirichlet_reads.py"


rule get_grey_matter_reads_se:
    input:
        fastq=config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.fastq.gz",
        dirichlet_matrix=config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
    output:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet.fastq.gz",
    log:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_SE.log",
    benchmark:
        repeat("benchmarks/get_{sample}_Grey_Matter.benchmark.txt", 1)
    message:
        "Preparing fastq files with all the reads that got assigned to the Grey Matter "
        "for sample {wildcards.sample} {MESSAGE_SUFFIX}"
    conda:
        "../envs/grey_matter_reads.yaml"
    script:
        "../scripts/get_grey_matter_reads_se.py"


rule get_grey_matter_reads_pe:
    input:
        fastq_r1=config["analysis_output_dir"] + "/fastq/PE/{sample}_R1_mapq.fastq.gz",
        fastq_r2=config["analysis_output_dir"] + "/fastq/PE/{sample}_R2_mapq.fastq.gz",
        dirichlet_matrix=config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
    output:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_R1.fastq.gz",
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_R2.fastq.gz",
    log:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_SE.log",
    benchmark:
        repeat("benchmarks/get_{sample}_Grey_Matter.benchmark.txt", 1)
    message:
        "Preparing fastq files with all the reads that got assigned to the Grey Matter "
        "for sample {wildcards.sample} {MESSAGE_SUFFIX}"
    conda:
        "../envs/grey_matter_reads.yaml"
    script:
        "../scripts/get_grey_matter_reads_pe.py"


rule get_dark_matter_reads_se:
    input:
        fastq=get_inputs_for_bowtie_r1,
        dirichlet_matrix=config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
    output:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet.fastq.gz",
    log:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet_SE.log",
    benchmark:
        repeat("benchmarks/get_{sample}_Dark_Matter.benchmark.txt", 1)
    message:
        "Preparing fastq files with all the reads that got assigned to the Dark Matter "
        "for sample {wildcards.sample} {MESSAGE_SUFFIX}"
    conda:
        "../envs/grey_matter_reads.yaml"
    script:
        "../scripts/get_dark_matter_reads_se.py"


rule get_dark_matter_reads_pe:
    input:
        fastq_r1=get_inputs_for_bowtie_r1,
        fastq_r2=get_inputs_for_bowtie_r2,
        dirichlet_matrix=config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
    output:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet_R1.fastq.gz",
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet_R2.fastq.gz",
    log:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet_SE.log",
    benchmark:
        repeat("benchmarks/get_{sample}_Dark_Matter.benchmark.txt", 1)
    message:
        "Preparing fastq files with all the reads that got assigned to the Dark Matter "
        "for sample {wildcards.sample} {MESSAGE_SUFFIX}"
    conda:
        "../envs/grey_matter_reads.yaml"
    script:
        "../scripts/get_dark_matter_reads_pe.py"


# noinspection PyUnresolvedReferences
def get_dirichlet_bams(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """

    sequences = get_total_paths(
        wildcards,
        checkpoints,
        config["with_entrez_query"],
        config["refseq_rep"],
        config["with_custom_sequences"],
        config["with_custom_accessions"],
        config["genera"],
        config["accessions"],
        config["sequences"],
    )

    inputs = []

    reads = ""
    if config["PE_MODERN"]:
        reads = "PE"
    elif config["PE_ANCIENT"] or config["SE"]:
        reads = "SE"

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["GBSeq_accession-version"],
        )

        inputs.append(
            config["analysis_output_dir"]
            + "/dirichlet_reads/{sample}/{orgname}/{orgname}_{accession}_dirichlet_{reads}.bam".format(
                sample=wildcards.sample, orgname=orgname, accession=accession, reads=reads,
            )
        )

    return inputs


def get_grey_matter_reads():

    if config["PE_MODERN"]:
        return [
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_R1.fastq.gz",
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_R2.fastq.gz",
        ]

    elif config["PE_ANCIENT"] or config["SE"]:
        return config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet.fastq.gz"


def get_dark_matter_reads():

    if config["PE_MODERN"]:
        return [
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet_R1.fastq.gz",
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet_R2.fastq.gz",
        ]

    elif config["PE_ANCIENT"] or config["SE"]:
        return config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet.fastq.gz"


# get dark matter, get grey matter

# input all


rule all_dirichlet:
    input:
        get_dirichlet_bams,
        get_grey_matter_reads(),
        get_dark_matter_reads(),
    output:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}_dirichlet_reads.done",
    benchmark:
        repeat("benchmarks/all_dirichlet_reads_{sample}.benchmark.txt", 1)
    message:
        "All the dirichlet assigned reads have been put in bam files for {wildcards.sample}."
    shell:
        "touch {output}"
