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
          "../scripts/count_accession_ts_tv.py"  # TODO keep to the naming scheme of {rulefile}_{rulename} as uses elsewhere


rule initial_ts_tv:
    input:
        # TODO replace with an input function that returns all the files to concatenate
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}.csv"
    output:
        # TODO don't use fake output names!
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}_aggregating.done"
    log:
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}_aggregating.log"
    params:
        query="{query}",  # TODO these are unnecesary, as you can access them with {wildcards.query}
        sample="{sample}"
    shell:
        "cat {input} 1>> {params.query}/ts_tv_counts/{params.sample}/allTsTvCounts.csv 2> {log}; touch {output}"


rule make_matrices_equal:  # TODO needs more literal name... I can't tell what this does from the current name
    input:
        "{query}/ts_tv_counts/{sample}/allTsTvCounts.csv",  # TODO don't use CamelCase in filenames (be consistent!)
        "{query}/fastq/{sample}_mapq.readlen",
        "{query}/entrez/{query}-selected-seqs.tsv"
    output:
        "{query}/probabilities/{sample}/TsTvMatrix.csv",
        "{query}/probabilities/{sample}/Probability_Model_Params.csv"  # TODO make this a .json file
    log:
        "{query}/probabilities/{sample}/TsTvMatrix.log"
    script:
        "../scripts/make_matrices_equal.py"  # TODO keep to the naming scheme


rule calculate_probabilities:
    input:
        "{query}/probabilities/{sample}/TsTvMatrix.csv",
        "{query}/probabilities/{sample}/Probability_Model_Params.csv",
        "fastq/{sample}.size"
    output:
        "{query}/probabilities/{sample}/Posterior_Probabilities.csv"  # TODO keep to the naming scheme
    log:
        "{query}/probabilities/{sample}/Posterior_Probabilities.log"
    params:
        submatrices=False  # TODO what is this?
    script:
          "../scripts/calculate_taxa_probabilities.py"  # TODO keep to the naming scheme


rule coverage_t_test:
    input:
        "{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam",
        "{query}/entrez/{query}-selected-seqs.tsv",  # TODO drop the dependencies on both these tsv files
        "{query}/entrez/{query}-nuccore.tsv"
    output:
        "{query}/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}.txt"
    log:
        "{query}/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}.log"
    script:
        "../scripts/coverage_t_test.py"


rule cat_pvalues:
    # TODO you must have an input: for this rule which contains all the outputs from coverage_t_test
    output:
        "{query}/probabilities/{sample}/t_test_pvalues.txt",
    log:
        "{query}/probabilities/{sample}/t_test_pvalues_agg.log"
    params:
        query="{query}",  # TODO these are unnecessary
        sample="{sample}"
    shell:
         # TODO it doesn't find all the files to concat because is has no inputs!!!
         #      so the scheduler will try to run it before any files exists as it has no depedencies
         "cat {params.query}/probabilities/{params.sample}/*_t_test_pvalue_*.txt 1>> {output} 2>> {log}"


rule calculate_dirichlet_abundances:
    input:
        "{query}/probabilities/{sample}/TsTvMatrix.csv",
        "{query}/probabilities/{sample}/t_test_pvalues.txt",
        "fastq/{sample}.size"
    output:
        "{query}/probabilities/{sample}/{sample}_posterior_abundance.tsv"
    log:
        "{query}/probabilities/{sample}/{sample}_posterior_abundance.log"
    script:
        "../scripts/calculate_dirichlet_abundances.py"

