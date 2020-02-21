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
