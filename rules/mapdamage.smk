#!/usr/bin/env python
# -*- coding: utf-8 -*-

WITH_REFSEQ_REP = config["WITH_REFSEQ_REP"]
WITH_ENTREZ_QUERY = config["WITH_ENTREZ_QUERY"]
SPECIFIC_GENUS = config["SPECIFIC_GENUS"]
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

        refseq_genomes = pd.read_csv(refseq_rep_prok.output[0], sep="\t")
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

    inputs = []

    if SPECIFIC_GENUS:
        sequences = sequences[
            sequences["species"].str.contains("|".join(SPECIFIC_GENUS))
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


# raise RuntimeError('PE data, mapDamage cannot run with that input format. Either collapse the reads, '
#                    'use SE data or do not include that rule.')


rule all_mapdamage:
    input:
        get_mapdamage_out_dir_paths,
    output:
        "{query}/mapdamage/{sample}_mapdamage.done",
    benchmark:
        repeat("benchmarks/all_alignments_{query}_{sample}.benchmark.txt", 1)
    shell:
        "touch {output}"
