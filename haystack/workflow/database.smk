#!/usr/bin/env python
# -*- coding: utf-8 -*-


include: "rules/entrez.smk"
include: "rules/bowtie_index.smk"
include: "rules/bowtie_meta_index.smk"
include: "rules/entrez_custom_db.smk"
include: "rules/refseq.smk"


wildcard_constraints:
    query="[\w]+",
    sample="[\w]+",
    orgname="[\w.-]+",
    accession="[\w.-]+",
    chunk="\d+",
