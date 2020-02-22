#!/usr/bin/env python
# -*- coding: utf-8 -*-


##### Target rules #####

rule entrez_nuccore_query:
    output:
        "entrez/{query}/{query}-nuccore.tsv"
    log:
        "entrez/{query}/{query}-nuccore.log"
    script:
        "../scripts/entrez_nuccore_query.py"

rule entrez_taxa_query:
    input:
        "entrez/{query}/{query}-nuccore.tsv"
    output:
        "entrez/{query}/{query}-taxa.tsv"
    log:
        "entrez/{query}/{query}-taxa.log"
    script:
        "../scripts/entrez_taxonomy_query.py"

rule entrez_pick_sequences:
    input:
         "entrez/{query}/{query}-nuccore.tsv",
         "entrez/{query}/{query}-taxa.tsv"
    output:
         "entrez/{query}/{query}-selected-seqs.tsv"
    log:
         "entrez/{query}/{query}-selected-seqs.log"
    script:
        "../scripts/entrez_select_seqs.py"

# TODO Can I not have just a folder as an output, do I need a wildcard ?

rule entrez_dowload_sequenes:
    input:
         "entrez/{query}/{query}-selected-seqs.tsv"
    log:
         "pipeline_files/{query}-selected-seqs.log"
    output:
         directory("pipeline_files/{query}/refseq/")
    script:
         "../scripts/entrez_download_seqs.py"

rule make_bowtie_general_profile:
    input:
         directory("pipeline_files/{query}/refseq/")
    log:
         "pipeline_files/{query}-bowtie2-seqs.log"
    output:
          "pipeline_files/{query}/reference_genome/{query}.fasta"
    script:
          "../scripts/entrez_bowtie_multifasta.py"


