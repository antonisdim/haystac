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
    threads: 8
    conda:
        "../envs/samtools.yaml"
    shell:
        "path=$(awk -F'\t' '$1 == \"{wildcards.orgname}\" {{print $3}}' {input}); "
        'type=$(htsfile "$path"); '
        'if [[ "$type" == *"gzip-compressed"* ]]; then'
        '   bgzip --decompress --stdout --threads {threads} "$path" | bgzip --stdout --threads {threads} 1> {output}; '
        'elif [[ "$type" == *"BGZF-compressed"* ]]; then'
        '   cp "$path" {output}; '
        "else "
        '   bgzip --stdout --threads 8 "$path" 1> {output}; '
        "fi"
