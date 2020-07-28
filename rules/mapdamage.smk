#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


import pandas as pd

from scripts.rip_utilities import get_total_paths, normalise_name

##### Target rules #####


rule dedup_merged_mapdamage:
    input:
        bam=(
            config["analysis_output_dir"]
            + "/sigma/{sample}/{reads}/{orgname}/{orgname}_{accession}.bam"
        ),
    log:
        config[
            "analysis_output_dir"
        ] + "/mapdamage/{sample}/rmdup_bam/{reads}/{orgname}/{orgname}_{accession}_rmdup.log",
    output:
        config[
            "analysis_output_dir"
        ] + "/mapdamage/{sample}/rmdup_bam/{reads}/{orgname}/{orgname}_{accession}_rmdup.bam",
    benchmark:
        repeat(
            "benchmarks/dedup_merged_{sample}_{reads}_{orgname}_{accession}.benchmark.txt", 1,
        )
    params:
        output=(
            config["analysis_output_dir"] + "/mapdamage/{sample}/rmdup_bam/{reads}/{orgname}/"
        ),
    message:
        "Removing duplicate reads that were aligned against genome {wildcards.accession} "
        "for taxon {wildcards.orgname}, for sample {wildcards.sample}. The unique aligned reads can be found "
        "in {output} and the log file can be found in {log}."
    shell:
        "dedup --merged --input {input.bam} --output {params.output} &> {log}"


rule run_mapdamage:
    input:
        bam=(
            config["analysis_output_dir"]
            + "/mapdamage/{sample}/rmdup_bam/{reads}/{orgname}/{orgname}_{accession}_rmdup.bam"
        ),
        ref_genome=config["genome_cache_folder"] + "/{orgname}/{accession}.fasta.gz",
    log:
        config["analysis_output_dir"] + "/mapdamage/{sample}/{reads}/{orgname}_{accession}.log",
    output:
        directory(
            config["analysis_output_dir"] + "/mapdamage/{sample}/{reads}/{orgname}-{accession}"
        ),
    message:
        "Performing a mapDamage analysis on unique aligned reads against genome {wildcards.accession} "
        "for taxon {wildcards.orgname}, for sample {wildcards.sample}. The output can be found in {output}, "
        "and its log file can be found in {log}."
    shell:
        "mapDamage -i {input.bam} -r {input.ref_genome} -d {output} --merge-libraries 2> {log}"


# noinspection PyUnresolvedReferences
def get_mapdamage_out_dir_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """

    sequences = get_total_paths(
        wildcards,
        checkpoints,
        config["with_entrez_query"],
        config["refseq_rep"],
        config["with_custom_sequences"],
        config["with_custom_accessions"],
        config["genera"],
    )

    inputs = []

    reads = ""
    if config["PE_MODERN"]:
        reads = "PE"
    elif config["PE_ANCIENT"] or config["SE"]:
        reads = "SE"

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["GBSeq_accession-version"],
        )

        inputs.append(
            config["analysis_output_dir"]
            + "/mapdamage/{sample}/{reads}/{orgname}-{accession}".format(
                sample=wildcards.sample,
                orgname=orgname,
                accession=accession,
                reads=reads,
            )
        )

    if config["PE_ANCIENT"] or config["SE"]:
        return inputs
    else:
        print(
            "WARNING: dedup is treating PE uncollapsed reads as SE reads. "
            "Removing PCR duplicates might not have been done correctly."
        )
        print(
            "WARNING: mapDamage has not been optimised to analyse paired end alignment data."
        )
        return inputs


rule all_mapdamage:
    input:
        get_mapdamage_out_dir_paths,
    output:
        config["analysis_output_dir"] + "/mapdamage/{sample}_mapdamage.done",
    benchmark:
        repeat("benchmarks/all_alignments_{sample}.benchmark.txt", 1)
    message:
        "MapDamage analysis have been performed for all the taxa in our database for sample {wildcards.sample}."
    shell:
        "touch {output}"
