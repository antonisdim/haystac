#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

WITH_REFSEQ_REP = config["WITH_REFSEQ_REP"]  # TODO get rid of all this redundancy
WITH_ENTREZ_QUERY = config["WITH_ENTREZ_QUERY"]
WITH_CUSTOM_SEQUENCES = config["WITH_CUSTOM_SEQUENCES"]
WITH_CUSTOM_ACCESSIONS = config["WITH_CUSTOM_ACCESSIONS"]
SPECIFIC_GENERA = config["SPECIFIC_GENERA"]
PE_ANCIENT = config["PE_ANCIENT"]
PE_MODERN = config["PE_MODERN"]
SE = config["SE"]

import pandas as pd

##### Target rules #####


rule dedup_merged_mapdamage:
    input:
        bam="{query}/sigma/{sample}/{reads}/{orgname}/{orgname}_{accession}.bam",
    log:
        "{query}/mapdamage/{sample}/rmdup_bam/{reads}/{orgname}/{orgname}_{accession}_rmdup.log",
    output:
        "{query}/mapdamage/{sample}/rmdup_bam/{reads}/{orgname}/{orgname}_{accession}_rmdup.bam",
    benchmark:
        repeat(
            "benchmarks/dedup_merged_{query}_{sample}_{reads}_{orgname}_{accession}.benchmark.txt",
            1,
        )
    params:
        output="{query}/mapdamage/{sample}/rmdup_bam/{reads}/{orgname}/",
    message:
        "Removing duplicate reads that were aligned against genome {wildcards.accession} "
        "for taxon {wildcards.orgname}, for sample {wildcards.sample}. The unique aligned reads can be found "
        "in {output} and the log file can be found in {log}."
    shell:
        "dedup --merged --input {input.bam} --output {params.output} &> {log}"


rule run_mapdamage:
    input:
        bam="{query}/mapdamage/{sample}/rmdup_bam/{reads}/{orgname}/{orgname}_{accession}_rmdup.bam",
        ref_genome="database/{orgname}/{accession}.fasta.gz",
    log:
        "{query}/mapdamage/{sample}/{reads}/{orgname}_{accession}.log",
    output:
        directory("{query}/mapdamage/{sample}/{reads}/{orgname}-{accession}"),
    message:
        "Performing a mapDamage analysis on unique aligned reads against genome {wildcards.accession} "
        "for taxon {wildcards.orgname}, for sample {wildcards.sample}. The output can be found in {output}, "
        "and its log file can be found in {log}."
    shell:
        "mapDamage -i {input.bam} -r {input.ref_genome} -d {output}"


# noinspection PyUnresolvedReferences
def get_mapdamage_out_dir_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """

    sequences = pd.DataFrame()

    if WITH_ENTREZ_QUERY:
        pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
        sequences = pd.read_csv(pick_sequences.output[0], sep="\t")

        if len(sequences) == 0:
            raise RuntimeError("The entrez pick sequences file is empty.")

    if WITH_REFSEQ_REP:
        refseq_rep_prok = checkpoints.entrez_refseq_accessions.get(
            query=wildcards.query
        )

        refseq_genomes = pd.read_csv(
            refseq_rep_prok.output[0], sep="\t"
        )  # TODO use output names, not indicies
        genbank_genomes = pd.read_csv(refseq_rep_prok.output[1], sep="\t")
        assemblies = pd.read_csv(refseq_rep_prok.output[2], sep="\t")
        refseq_plasmids = pd.read_csv(refseq_rep_prok.output[3], sep="\t")
        genbank_plasmids = pd.read_csv(refseq_rep_prok.output[4], sep="\t")

        invalid_assemblies = checkpoints.entrez_invalid_assemblies.get(
            query=wildcards.query
        )
        invalid_assembly_sequences = pd.read_csv(invalid_assemblies.output[0], sep="\t")

        assemblies = assemblies[
            ~assemblies["GBSeq_accession-version"].isin(
                invalid_assembly_sequences["GBSeq_accession-version"]
            )
        ]

        # TODO get rid of the redundancy!! we shouldn't be duplicating code
        if WITH_ENTREZ_QUERY:
            sequences = pd.concat(
                [
                    sequences,
                    refseq_genomes,
                    genbank_genomes,
                    assemblies,
                    refseq_plasmids,
                    genbank_plasmids,
                ]
            )
        else:
            sequences = pd.concat(
                [
                    refseq_genomes,
                    genbank_genomes,
                    assemblies,
                    refseq_plasmids,
                    genbank_plasmids,
                ]
            )

    if WITH_CUSTOM_SEQUENCES:
        custom_fasta_paths = pd.read_csv(
            config["custom_seq_file"],
            sep="\t",
            header=None,
            names=["species", "GBSeq_accession-version", "path"],
        )

        custom_seqs = custom_fasta_paths[["species", "GBSeq_accession-version"]]

        sequences = sequences.append(custom_seqs)

    if WITH_CUSTOM_ACCESSIONS:
        custom_accessions = pd.read_csv(
            config["custom_acc_file"],
            sep="\t",
            header=None,
            names=["species", "GBSeq_accession-version"],
        )

        sequences = sequences.append(custom_accessions)

    inputs = []

    if SPECIFIC_GENERA:
        sequences = sequences[
            sequences["species"].str.contains("|".join(SPECIFIC_GENERA))
        ]

    reads = ""
    if PE_MODERN:
        reads = "PE"
    elif PE_ANCIENT or SE:
        reads = "SE"

    for key, seq in sequences.iterrows():
        orgname, accession = (
            seq["species"].replace(" ", "_").replace("[", "").replace("]", ""),
            seq["GBSeq_accession-version"],
        )

        inputs.append(
            "{query}/mapdamage/{sample}/{reads}/{orgname}-{accession}".format(
                query=wildcards.query,
                sample=wildcards.sample,
                orgname=orgname,
                accession=accession,
                reads=reads,
            )
        )

    if PE_ANCIENT or SE:
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


# TODO delete commented out code
# raise RuntimeError('PE data, mapDamage cannot run with that input format. Either collapse the reads, '
#                    'use SE data or do not include that rule.')


# TODO delete ruls that are not used by the pipeline
rule all_mapdamage:
    input:
        get_mapdamage_out_dir_paths,
    output:
        "{query}/mapdamage/{sample}_mapdamage.done",
    benchmark:
        repeat("benchmarks/all_alignments_{query}_{sample}.benchmark.txt", 1)
    message:
        "MapDamage analysis have been performed for all the taxa in our database for sample {wildcards.sample}."
    shell:
        "touch {output}"
