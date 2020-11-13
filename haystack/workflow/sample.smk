#!/usr/bin/env python
# -*- coding: utf-8 -*-


include: "rules/sample.smk"
include: "rules/adapterremoval.smk"
include: "rules/sra.smk"


wildcard_constraints:
    query="[\w]+",
    sample="[\w]+",
    orgname="[\w.-]+",
    accession="[\w.-]+",
    chunk="\d+",
