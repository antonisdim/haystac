#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

if config["module"] == "sample":

    include: "rules/sample.smk"

    include: "rules/adapterremoval.smk"

    include: "rules/sra.smk"


elif config["module"] == "database":

    include: "rules/entrez.smk"

    include: "rules/bowtie_index.smk"

    include: "rules/bowtie_meta_index.smk"

    include: "rules/entrez_custom_db.smk"

    include: "rules/refseq.smk"


elif config["module"] == "analyse":

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
