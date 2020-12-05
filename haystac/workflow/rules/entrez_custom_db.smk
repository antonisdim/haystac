#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


rule entrez_custom_sequences:
    input:
        config["sequences"],
    log:
        config["cache"] + "/ncbi/{orgname}/custom_seq-{accession}.log",
    output:
        config["cache"] + "/ncbi/{orgname}/custom_seq-{accession}.fasta.gz",
    message:
        "Adding the user provided fasta sequence {wildcards.accession} for taxon {wildcards.orgname} to the database."
    shell:
        "(grep {wildcards.orgname} {input} | cut -f 3 | xargs -I '{{}}' bash -c '[[ '{{}}' == *.gz ]] && "
        "bgzip -cd '{{}}' 1> " + config["cache"] + "/ncbi/{wildcards.orgname}/temp-{wildcards.accession}.fasta; "
        "bgzip -c " + config["cache"] + "/ncbi/{wildcards.orgname}/temp-{wildcards.accession}.fasta 1> {output} || "
        "bgzip -c '{{}}' 1> {output}') 2> {log}"
