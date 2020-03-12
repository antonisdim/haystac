#!/usr/bin/env python
# -*- coding: utf-8 -*-

##### Target rules #####


rule parse_bams:
    input:
        "{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam"
    output:
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}.csv"
    log:
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}.log"
    params:
        taxon="{orgname}"
    script:
          "../scripts/parse_bams.py"



rule initial_ts_tv:
    input:
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}.csv"
    output:
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}_aggregating.done"
    log:
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}_aggregating.log"
    params:
        query="{query}",
        sample="{sample}"
    shell:
        "cat {input} 1>> {params.query}/ts_tv_counts/{params.sample}/allTsTvCounts.csv 2> {log}; touch {output}"



rule make_matrices_equal:
    input:
        "{query}/ts_tv_counts/{sample}/allTsTvCounts.csv",
        "{query}/fastq/{sample}_mapq.readlen",
        "{query}/entrez/{query}-selected-seqs.tsv"
    output:
        "{query}/probabilities/{sample}/TsTvMatrix.csv",
        "{query}/probabilities/{sample}/Probability_Model_Params.csv"
    log:
        "{query}/probabilities/{sample}/TsTvMatrix.log"
    script:
        "../scripts/make_matrices_equal.py"



rule calculate_probabilities:
    input:
        "{query}/probabilities/{sample}/TsTvMatrix.csv",
        "{query}/probabilities/{sample}/Probability_Model_Params.csv",
        "fastq/{sample}.size"
    output:
        "{query}/probabilities/{sample}/Posterior_Probabilities.csv"
    log:
        "{query}/probabilities/{sample}/Posterior_Probabilities.log"
    params:
        sample_name="{sample}",
        submatrices=False
    script:
          "../scripts/calculate_taxa_probabilities.py"



rule coverage_t_test:
    input:
        "{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam",
        "{query}/entrez/{query}-selected-seqs.tsv",
        "{query}/entrez/{query}-nuccore.tsv"
    output:
        "{query}/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}.txt"
    log:
        "{query}/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}.log"
    params:
        taxon="{orgname}"
    script:
        "../scripts/coverage_t_test.py"



rule cat_pvalues:
    output:
        "{query}/probabilities/{sample}/t_test_pvalues.txt",
    log:
        "{query}/probabilities/{sample}/t_test_pvalues_agg.log"
    params:
        query="{query}",
        sample="{sample}"
    shell:
         # todo it is dodgy as somtimes it doesn't find all the files to concat
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

