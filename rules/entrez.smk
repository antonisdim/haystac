#!/usr/bin/env python
# -*- coding: utf-8 -*-



##### Target rules #####

rule entrez_nuccore_query:
    output:
        "entrez/{query}/{query}-nuccore.tsv"
    log:
        "entrez/{query}/{query}-nuccore.log"
    script:
        "../scripts/get_nuccore_query.py"
