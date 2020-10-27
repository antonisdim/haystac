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
    orgname="[^/]+", # TODO check for consistency with other parts of the code
    accession="[^/]+", # TODO check for consistency with other parts of the code
    chunk="\d+",
