#!/usr/bin/env python
# -*- coding: utf-8 -*-


include: "rules/entrez.smk"
include: "rules/bowtie.smk"
include: "rules/bowtie_meta.smk"
include: "rules/metagenomics.smk"
include: "rules/entrez_custom_db.smk"
include: "rules/refseq.smk"
include: "rules/adapterremoval.smk"
include: "rules/mapdamage.smk"
include: "rules/sra.smk"
include: "rules/split_reads.smk"


wildcard_constraints:
    query="[\w]+",
    sample="[\w]+",
    orgname="[\w.-]+",
    accession="[\w.-]+",
    chunk="\d+",
