#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from haystack.workflow.scripts.utilities import get_total_paths, normalise_name

MESSAGE_SUFFIX = "(output: {output} and log: {log})" if config["debug"] else ""


def get_bams_for_ts_tv_count(wildcards):
    if config["PE_MODERN"]:
        return config["analysis_output_dir"] + "/alignments/{sample}/PE/{orgname}/{orgname}_{accession}.bam".format(
            sample=wildcards.sample, orgname=wildcards.orgname, accession=wildcards.accession,
        )
    elif config["PE_ANCIENT"] or config["SE"]:
        return config["analysis_output_dir"] + "/alignments/{sample}/SE/{orgname}/{orgname}_{accession}.bam".format(
            sample=wildcards.sample, orgname=wildcards.orgname, accession=wildcards.accession,
        )


rule count_accession_ts_tv:
    input:
        get_bams_for_ts_tv_count,
    output:
        config["analysis_output_dir"] + "/ts_tv_counts/{sample}/{orgname}_count_{accession}.csv",
    log:
        config["analysis_output_dir"] + "/ts_tv_counts/{sample}/{orgname}_count_{accession}.log",
    benchmark:
        repeat(
            "benchmarks/count_accession_ts_tv_{sample}_{orgname}_{accession}.benchmark.txt", 1,
        )
    params:
        pairs=config["PE_MODERN"],
    message:
        "Counting the number of transitions and transversions per read for taxon {wildcards.orgname} {MESSAGE_SUFFIX}"
    script:
        "../scripts/count_accession_ts_tv.py"


def get_ts_tv_count_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """
    sequences = get_total_paths(
        wildcards,
        checkpoints,
        config["query"],
        config["refseq_rep"],
        config["sequences"],
        config["accessions"],
        config["genera"],
    )

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["GBSeq_accession-version"],
        )

        inputs.append(
            config["analysis_output_dir"]
            + "/ts_tv_counts/{sample}/{orgname}_count_{accession}.csv".format(
                sample=wildcards.sample, orgname=orgname, accession=accession,
            )
        )

    return inputs


rule initial_ts_tv:
    input:
        get_ts_tv_count_paths,
    output:
        config["analysis_output_dir"] + "/ts_tv_counts/{sample}/all_ts_tv_counts.csv",
    log:
        config["analysis_output_dir"] + "/ts_tv_counts/{sample}/all_ts_tv_counts.log",
    benchmark:
        repeat("benchmarks/initial_ts_tv_{sample}.benchmark.txt", 1)
    message:
        "Concatenating all the Ts and Tv count files for sample {wildcards.sample} {MESSAGE_SUFFIX}"
    script:
        "../scripts/concat_files.py"


def get_right_readlen(wildcards):
    if config["PE_MODERN"]:
        return config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq_pair.readlen".format(sample=wildcards.sample)
    else:
        return config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.readlen".format(sample=wildcards.sample)


rule calculate_likelihoods:
    input:
        config["analysis_output_dir"] + "/ts_tv_counts/{sample}/all_ts_tv_counts.csv",
        get_right_readlen,
        get_ts_tv_count_paths,
    output:
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_probability_model_params.json",
    log:
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.log",
    benchmark:
        repeat("benchmarks/calculate_likelihoods_{sample}.benchmark.txt", 1)
    message:
        "Calculating the likelihoods and performing the Dirichlet assignment of the reads in sample "
        "{wildcards.sample} to the taxa in our database {MESSAGE_SUFFIX}"
    script:
        "../scripts/calculate_likelihoods.py"


rule calculate_taxa_probabilities:
    input:
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_probability_model_params.json",
        config["sample_output_dir"] + "/fastq_inputs/meta/{sample}.size",
    output:
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_posterior_probabilities.tsv",
    log:
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_posterior_probabilities.log",
    benchmark:
        repeat("benchmarks/calculate_taxa_probabilities_{sample}.benchmark.txt", 1)
    params:
        submatrices=False,
    message:
        "Calculating the taxonomic assignment posterior probabilities for sample {wildcards.sample} {MESSAGE_SUFFIX}"
    script:
        "../scripts/calculate_taxa_probabilities.py"


rule coverage_t_test:
    input:
        config["analysis_output_dir"] + "/alignments/{sample}/{reads}/{orgname}/{orgname}_{accession}.bam",
        config["cache"] + "/{orgname}/{accession}.fasta.gz.fai",
    output:
        config["analysis_output_dir"] + "/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}_{reads}.txt",
    log:
        config["analysis_output_dir"] + "/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}_{reads}.log",
    benchmark:
        repeat(
            "benchmarks/coverage_t_test_{sample}_{orgname}_{accession}_{reads}.benchmark.txt", 1,
        )
    message:
        "Performing a T-Test to assess if reads from sample {wildcards.sample} represent "
        "a random genome sample of taxon {wildcards.orgname} {MESSAGE_SUFFIX}"
    script:
        "../scripts/coverage_t_test.py"


def get_t_test_values_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """

    sequences = get_total_paths(
        wildcards,
        checkpoints,
        config["query"],
        config["refseq_rep"],
        config["sequences"],
        config["accessions"],
        config["genera"],
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
            + "/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}_{reads}.txt".format(
                sample=wildcards.sample, orgname=orgname, accession=accession, reads=reads,
            )
        )

    return inputs


rule cat_pvalues:
    input:
        get_t_test_values_paths,
    output:
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_t_test_pvalues.txt",
    log:
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_t_test_pvalues.log",
    benchmark:
        repeat("benchmarks/cat_pvalues_{sample}.benchmark.txt", 1)
    message:
        "Concatenating all the T-Test p-value outputs for sample {wildcards.sample} {MESSAGE_SUFFIX}"
    script:
        "../scripts/concat_files.py"


rule calculate_dirichlet_abundances:
    input:
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_t_test_pvalues.txt",
        config["sample_output_dir"] + "/fastq_inputs/meta/{sample}.size",
    output:
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_posterior_abundance.tsv",
    log:
        config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_posterior_abundance.log",
    benchmark:
        repeat("benchmarks/calculate_dirichlet_abundances_{sample}.benchmark.txt", 1)
    message:
        "Calculating the mean posterior abundance for sample {wildcards.sample} {MESSAGE_SUFFIX}"
    script:
        "../scripts/calculate_dirichlet_abundances.py"
