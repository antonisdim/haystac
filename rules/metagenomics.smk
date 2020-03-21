#!/usr/bin/env python
# -*- coding: utf-8 -*-


##### Target rules #####


rule count_fastq_length:
    input:
         fastq=lambda wildcards: config['samples'][wildcards.sample]
    log:
         "fastq/{sample}.log"
    output:
         "fastq/{sample}.size"
    shell:
          "seqtk seq -A {input.fastq} | grep -v '^>' | wc -l 1> {output} 2> {log}"



rule count_accession_ts_tv:
    input:
        "{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam"
    output:
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}.csv"
    log:
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}.log"
    script:
          "../scripts/count_accession_ts_tv.py"



# noinspection PyUnresolvedReferences
def get_ts_tv_count_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """
    pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
    sequences = pd.read_csv(pick_sequences.output[0], sep='\t')

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = seq['species'].replace(" ", "_"), seq['GBSeq_accession-version']

        inputs.append('{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}.csv'.
                      format(query=wildcards.query, sample=wildcards.sample, orgname=orgname, accession=accession))

    return inputs



rule initial_ts_tv:
    input:
       get_ts_tv_count_paths
    output:
        "{query}/ts_tv_counts/{sample}/all_ts_tv_counts.csv"
    log:
        "{query}/ts_tv_counts/{sample}/all_ts_tv_counts.log"
    shell:
        "cat {input} 1> {output} 2> {log}"



rule calculate_likelihoods:
    input:
        "{query}/ts_tv_counts/{sample}/all_ts_tv_counts.csv",
        "{query}/fastq/{sample}_mapq.readlen",
        "{query}/entrez/{query}-selected-seqs.tsv"
    output:
        "{query}/probabilities/{sample}/likelihood_ts_tv_matrix.csv",
        "{query}/probabilities/{sample}/probability_model_params.json"
    log:
        "{query}/probabilities/{sample}/likelihood_ts_tv_matrix.log"
    script:
        "../scripts/calculate_likelihoods.py" #todo ask Evan to check if they are the same with the SQL commands



rule calculate_probabilities:
    input:
        "{query}/probabilities/{sample}/likelihood_ts_tv_matrix.csv",
        "{query}/probabilities/{sample}/probability_model_params.json",
        "fastq/{sample}.size"
    output:
        "{query}/probabilities/{sample}/posterior_probabilities.csv"
    log:
        "{query}/probabilities/{sample}/posterior_probabilities.log"
    params:
        submatrices=False  # TODO what is this? - I kinda explain it in the calculate_taxa_probabilities.py file
    script:
          "../scripts/calculate_taxa_probabilities.py"  # TODO keep to the naming scheme


rule coverage_t_test:
    input:
        "{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam",
        "{query}/entrez/{query}-selected-seqs.tsv",  # TODO drop the dependencies on both these tsv files - how ? I need these files
        "{query}/entrez/{query}-nuccore.tsv"
    output:
        "{query}/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}.txt"
    log:
        "{query}/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}.log"
    script:
        "../scripts/coverage_t_test.py"



# noinspection PyUnresolvedReferences
def get_t_test_values_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """
    pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
    sequences = pd.read_csv(pick_sequences.output[0], sep='\t')

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = seq['species'].replace(" ", "_"), seq['GBSeq_accession-version']

        inputs.append('{query}/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}.txt'.
                      format(query=wildcards.query, sample=wildcards.sample, orgname=orgname, accession=accession))

    return inputs



rule cat_pvalues:
    input:
        get_t_test_values_paths
    output:
        "{query}/probabilities/{sample}/t_test_pvalues.txt",
    log:
        "{query}/probabilities/{sample}/t_test_pvalues.log"
    shell:
         "cat {input} 1> {output} 2> {log}"



rule calculate_dirichlet_abundances:
    input:
        "{query}/probabilities/{sample}/likelihood_ts_tv_matrix.csv",
        "{query}/probabilities/{sample}/t_test_pvalues.txt",
        "fastq/{sample}.size"
    output:
        "{query}/probabilities/{sample}/{sample}_posterior_abundance.tsv"
    log:
        "{query}/probabilities/{sample}/{sample}_posterior_abundance.log"
    script:
        "../scripts/calculate_dirichlet_abundances.py"

