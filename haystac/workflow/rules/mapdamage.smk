#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from haystac.workflow.scripts.utilities import get_total_paths, PE, print_warning


rule dedup_merged_mapdamage:
    input:
        bam=(
            config["analysis_output_dir"]
            + "/dirichlet_reads/{sample}/{orgname}/{orgname}_{accession}_dirichlet_{reads}.bam"
        ),
    log:
        config[
            "analysis_output_dir"
        ] + "/rmdup_bam/{sample}/{reads}/{orgname}/{orgname}_{accession}_dirichlet_{reads}_rmdup.log",
    output:
        config[
            "analysis_output_dir"
        ] + "/rmdup_bam/{sample}/{reads}/{orgname}/{orgname}_{accession}_dirichlet_{reads}_rmdup.bam",
    benchmark:
        repeat(
            "benchmarks/dedup_merged_{sample}_{reads}_{orgname}_{accession}.benchmark.txt", 1,
        )
    params:
        output=config["analysis_output_dir"] + "/rmdup_bam/{sample}/{reads}/{orgname}/",
    message:
        "Removing duplicate reads that were aligned to taxon {wildcards.orgname}, for sample {wildcards.sample}."
    conda:
        "../envs/dedup.yaml"
    shell:
        "dedup --merged --input {input.bam} --output {params.output} &> {log}"


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
        directory(config["analysis_output_dir"] + "/mapdamage/{sample}/{reads}/{orgname}-{accession}"),
    message:
        "Performing a mapDamage analysis for taxon {wildcards.orgname}, for sample {wildcards.sample}."
    conda:
        "../envs/mapdamage.yaml"
    shell:
        "mapDamage -i {input.bam} -r {input.ref_genome} -d {output} --merge-libraries 2> {log}"


# noinspection PyUnresolvedReferences
def get_mapdamage_out_dir_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """

    sequences = get_total_paths(checkpoints, config)
    inputs = []

    for orgname, accession in sequences:
        inputs.append(
            config["analysis_output_dir"] + f"/mapdamage/{wildcards.sample}/{config['read_mode']}/{orgname}-{accession}"
        )

    if config["read_mode"] == PE:
        print_warning("dedup treats uncollapsed PE reads as SE. PCR deduplication might not have been done correctly.")
        print_warning("mapDamage has not been optimised to analyse paired end alignment data.")

    return inputs


rule all_mapdamage:
    input:
        get_mapdamage_out_dir_paths,
    output:
        config["analysis_output_dir"] + "/mapdamage/{sample}_mapdamage.done",
    benchmark:
        repeat("benchmarks/all_alignments_{sample}.benchmark.txt", 1)
    message:
        "A mapDamage analysis has been performed for all the taxa in the database for sample {wildcards.sample}."
    shell:
        "touch {output}"
