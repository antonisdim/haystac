#!/usr/bin/env python
# -*- coding: utf-8 -*-

configfile: "config.yaml"

##### Modules #####

include: "rules/entrez.smk"
include: "rules/bowtie.smk"
include: "rules/bowtie_meta.smk"
include: "rules/metagenomics.smk"
include: "rules/entrez_alt.smk"
include: "rules/entrez_build_prok_refseq_rep.smk"
include: "rules/data_preprocessing.smk"
include: "rules/mapdamage.smk"

##### Wildcards #####

wildcard_constraints:
    query="[\w]+",
    sample="[\w]+",
    # orgname="[\w]+", #todo these wildcards need to have numbers, letters and symbols like '_, .'
    orgname="[^/]+",
    accession="[^/]+",
    chunk='\d+',
    # sample_accession="^(.(?!(_R1|_R2)))*."
    ##### Target rules #####

from datetime import datetime

startTime = datetime.now()
rule all:
    input:
        # PE modern
        # "fastq_inputs/PE_mod/SRR1031215_R1_adRm.fastq.gz",  # test data_preprocessing.smk for merged ancient PE
        # "refseq_rep/entrez/refseq_rep-invalid-assemblies.tsv",
        # "refseq_rep/bowtie/refseq_rep_refseq_prok.fasta.gz", # test entrez_build_prok_refseq_rep.smk
        # "refseq_rep/bowtie/refseq_rep_entrez.fasta.gz", # test entrez.smk,
        # "refseq_rep/fastq/PE/SRR1031215_mapq_pair.readlen", # test bowtie.smk
        # "refseq_rep/sigma/SRR1031215_alignments.done",   # test bowtie_meta.smk,
        # "refseq_rep/probabilities/SRR1031215/SRR1031215_posterior_probabilities.csv", #test metagenomics.smk - probabilities
        # "refseq_rep/probabilities/SRR1031215/SRR1031215_posterior_abundance.tsv", #test metagenomics.smk - abundances
        # "refseq_rep/mapdamage/SRR1031215_mapdamage.done" # test mapdamage.smk
        # PE ancient
        # "fastq_inputs/PE_anc/SRR1031289_adRm.fastq.gz", # test data_preprocessing.smk for merged ancient PE
        # "refseq_rep/bowtie/refseq_rep_refseq_prok.fasta.gz", # test entrez_build_prok_refseq_rep.smk
        # "refseq_rep/bowtie/refseq_rep_entrez.fasta.gz", # test entrez.smk,
        # "refseq_rep/fastq/SE/SRR1031289_mapq.readlen", # test bowtie.smk
        # "refseq_rep/sigma/SRR1031289_alignments.done",   # test bowtie_meta.smk
        # "refseq_rep/probabilities/SRR1031289/SRR1031289_posterior_probabilities.csv", #test metagenomics.smk - probabilities
        # "refseq_rep/probabilities/SRR1031289/SRR1031289_posterior_abundance.tsv", #test metagenomics.smk - abundances
        # "refseq_rep/mapdamage/SRR1031289_mapdamage.done" # test mapdamage.smk
        # SE
        # "fastq_inputs/SE/SRR054920_adRm.fastq.gz",  # test data_preprocessing.smk for ancient SE
        # "refseq_rep/bowtie/refseq_rep_refseq_prok.fasta.gz",  # test entrez_build_prok_refseq_rep.smk
        # "refseq_rep/bowtie/refseq_rep_entrez.fasta.gz",  # test entrez.smk
        # "refseq_rep/fastq/SE/SRR054920_mapq.readlen",  # test bowtie.smk
        # "refseq_rep/sigma/SRR054920_alignments.done",  # test bowtie_meta.smk
        # "refseq_rep/probabilities/SRR054920/SRR054920_posterior_probabilities.csv",
        # #test metagenomics.smk - probabilities
        # "refseq_rep/probabilities/SRR054920/SRR054920_posterior_abundance.tsv",  #test metagenomics.smk - abundances
        # "refseq_rep/mapdamage/SRR054920_mapdamage.done" # test mapdamage.smk

        # "refseq_rep/bowtie/refseq_rep_refseq_prok.fasta.gz", # test entrez_build_prok_refseq_rep.smk
        # "refseq_rep/bowtie/refseq_rep_entrez.fasta.gz", # test entrez.smk
        # "refseq_rep/fastq/RISE00_mapq.readlen",  # test bowtie.smk
        # "refseq_rep/sigma/RISE00_alignments.done",   # test bowtie_meta.smk
        # "refseq_rep/probabilities/RISE00/RISE00_posterior_probabilities.csv", #test metagenomics.smk - probabilities
        # "refseq_rep/probabilities/RISE00/RISE00_posterior_abundance.tsv" #test metagenomics.smk - abundances
        # "yersinia_test/entrez_alt/sizes.txt",

        # "yersinia_test/bowtie/yersinia_test.fasta.gz",  # test entrez.smk
        # "yersinia_test/fastq/RISE00_mapq.readlen",      # test bowtie.smk
        # "yersinia_test/sigma/RISE00_alignments.done",   # test bowtie_meta.smk
        # "yersinia_test/probabilities/RISE00/RISE00_posterior_probabilities.csv", #test metagenomics.smk - probabilities
        # "yersinia_test/probabilities/RISE00/RISE00_posterior_abundance.tsv" #test metagenomics.smk - abundances

print(datetime.now() - startTime)
