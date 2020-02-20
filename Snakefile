#!/usr/bin/env python
# -*- coding: utf-8 -*-

from snakemake.utils import format

configfile: "config.yaml"

##### Modules #####

include: "rules/example.smk"

##### Wildcards #####

wildcard_constraints:
    reference="[\w.]+",
    n='\d+'


##### Target rules #####

rule all:
    input:
        "report/{reference}/file.txt"
