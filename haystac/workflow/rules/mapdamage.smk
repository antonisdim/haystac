#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from haystac.workflow.scripts.utilities import get_final_db_paths, PE, print_warning, print_error


rule dedup_merged_mapdamage:
    input:
        bam=(
            config["analysis_output_dir"]
            + "/dirichlet_reads/{sample}/{orgname}/{orgname}_{accession}_dirichlet_COLLAPSED.bam"
        ),
    log:
        config[
            "analysis_output_dir"
        ] + "/rmdup_bam/{sample}/COLLAPSED/{orgname}/{orgname}_{accession}_dirichlet_COLLAPSED_rmdup.log",
    output:
        config[
            "analysis_output_dir"
        ] + "/rmdup_bam/{sample}/COLLAPSED/{orgname}/{orgname}_{accession}_dirichlet_COLLAPSED_rmdup.bam",
    params:
        output=config["analysis_output_dir"] + "/rmdup_bam/{sample}/COLLAPSED/{orgname}/",
    message:
        "Removing duplicate reads that were aligned to taxon {wildcards.orgname}, for sample {wildcards.sample}."
    conda:
        "../envs/dedup.yaml"
    shell:
        "dedup --merged --input {input.bam} --output {params.output} &> {log}"


rule picard_single_mapdamage:
    input:
        bam=(
            config["analysis_output_dir"]
            + "/dirichlet_reads/{sample}/{orgname}/{orgname}_{accession}_dirichlet_{reads}.bam"
        ),
    log:
        config[
            "analysis_output_dir"
        ] + "/rmdup_bam/{sample}/SE/{orgname}/{orgname}_{accession}_dirichlet_{reads}_rmdup.log",
    output:
        config[
            "analysis_output_dir"
        ] + "/rmdup_bam/{sample}/SE/{orgname}/{orgname}_{accession}_dirichlet_{reads}_rmdup.bam",
    params:
        output=(
            config["analysis_output_dir"] + "/rmdup_bam/{sample}/SE/{orgname}/{orgname}_{accession}_dirichlet_{reads}"
        ),
    message:
        "Removing duplicate reads that were aligned to taxon {wildcards.orgname}, for sample {wildcards.sample}."
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates --INPUT {input.bam} --OUTPUT {params.output}_rmdup.bam "
        "--METRICS_FILE {params.output}_metrics.txt --REMOVE_DUPLICATES True &> {log}"


rule run_mapdamage:
    input:
        bam=(
            config["analysis_output_dir"]
            + "/rmdup_bam/{sample}/{reads}/{orgname}/{orgname}_{accession}_dirichlet_{reads}_rmdup.bam"
        ),
        ref_genome=config["cache"] + "/ncbi/{orgname}/{accession}.fasta.gz",
    log:
        config["analysis_output_dir"] + "/mapdamage/{sample}/{reads}/{orgname}_{accession}.log",
    output:
        directory(config["analysis_output_dir"] + "/mapdamage/{sample}/{reads}/{orgname}/{accession}"),
    message:
        "Performing a mapDamage analysis for taxon {wildcards.orgname}, for sample {wildcards.sample}."
    conda:
        "../envs/mapdamage.yaml"
    shell:
        "mapDamage -i {input.bam} -r {input.ref_genome} -d {output} 2> {log}"


# noinspection PyUnresolvedReferences
def get_mapdamage_out_dir_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """

    sequences = get_final_db_paths(checkpoints)
    inputs = []

    for orgname, accession in sequences:
        inputs.append(
            config["analysis_output_dir"] + f"/mapdamage/{wildcards.sample}/{config['read_mode']}/{orgname}/{accession}"
        )

    if config["read_mode"] == PE:
        print_error("mapDamage has not been optimised to analyse paired end alignment data.")

    return inputs


rule all_mapdamage:
    input:
        get_mapdamage_out_dir_paths,
    output:
        config["analysis_output_dir"] + "/mapdamage/{sample}_mapdamage.done",
    message:
        "A mapDamage analysis has been performed for all the taxa in the database for sample {wildcards.sample}."
    shell:
        "touch {output}"
