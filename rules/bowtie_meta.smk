#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pandas as pd

MIN_FRAG_LEN = 0  # I found them from the actual command that SIGMA runs
MAX_FRAG_LEN = 1000
SIGMA_MIN_SCORE_CONSTANT = (
    -6
)  # it's from SIGMA, wouldn't want this changed unless the user is SUPER knowledgeable
WITH_REFSEQ_REP = config["WITH_REFSEQ_REP"]
WITH_ENTREZ_QUERY = config["WITH_ENTREZ_QUERY"]

##### Target rules #####


rule index_database_entrez:
    input:
        "database/{orgname}/{accession}.fasta.gz",
    log:
        "database/{orgname}/{accession}_index.log",
    output:
        expand("database/{{orgname}}/{{accession}}.{n}.bt2l", n=[1, 2, 3, 4]),
        expand("database/{{orgname}}/{{accession}}.rev.{n}.bt2l", n=[1, 2]),
    benchmark:
        repeat("benchmarks/index_database_{orgname}_{accession}.benchmark.txt", 1)
    shell:
        "bowtie2-build --large-index {input} database/{wildcards.orgname}/{wildcards.accession} &> {log}"


def get_idx_entrez(wildcards):
    """
    Get all the index paths for the taxa in our database from the entrez query.
    """
    if not WITH_ENTREZ_QUERY:
        return []

    pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
    sequences = pd.read_csv(pick_sequences.output[0], sep="\t")

    if len(sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            seq["species"].replace(" ", "_").replace("[", "").replace("]", ""),
            seq["GBSeq_accession-version"],
        )
        inputs.append(
            "database/{orgname}/{accession}.1.bt2l".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


def get_idx_ref_gen(wildcards):
    """
    Get all the index paths for the taxa in our database from the refseq rep and genbank.
    """
    if not WITH_REFSEQ_REP:
        return []

    refseq_rep_prok = checkpoints.entrez_refseq_accessions.get(query=wildcards.query)

    refseq_genomes = pd.read_csv(refseq_rep_prok.output[0], sep="\t")
    genbank_genomes = pd.read_csv(refseq_rep_prok.output[1], sep="\t")
    refseq_plasmids = pd.read_csv(refseq_rep_prok.output[3], sep="\t")
    genbank_plasmids = pd.read_csv(refseq_rep_prok.output[4], sep="\t")

    sequences = pd.concat(
        [refseq_genomes, genbank_genomes, refseq_plasmids, genbank_plasmids]
    )

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            seq["species"].replace(" ", "_").replace("[", "").replace("]", ""),
            seq["GBSeq_accession-version"],
        )
        inputs.append(
            "database/{orgname}/{accession}.1.bt2l".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


def get_idx_assembly(wildcards):
    """
    Get all the individual bam file paths for the taxa in our database.
    """
    if not WITH_REFSEQ_REP:
        return []

    refseq_rep_prok = checkpoints.entrez_refseq_accessions.get(query=wildcards.query)

    assemblies = pd.read_csv(refseq_rep_prok.output[2], sep="\t")

    invalid_assemblies = checkpoints.entrez_invalid_assemblies.get(
        query=wildcards.query
    )
    invalid_assembly_sequences = pd.read_csv(invalid_assemblies.output[0], sep="\t")

    assemblies = assemblies[
        ~assemblies["GBSeq_accession-version"].isin(
            invalid_assembly_sequences["GBSeq_accession-version"]
        )
    ]

    sequences = assemblies

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            seq["species"].replace(" ", "_").replace("[", "").replace("]", ""),
            seq["GBSeq_accession-version"],
        )
        inputs.append(
            "database/{orgname}/{accession}.1.bt2l".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


rule idx_database:
    input:
        get_idx_entrez,
        get_idx_ref_gen,
        get_idx_assembly,
    log:
        "database/idx_database_{query}.log",
    output:
        "database/idx_database_{query}.done",
    shell:
        "touch {output}"


def get_min_score(wildcards, input):
    """Get the min score dor the edit distance of the alignment."""
    return (
        round(float(open(input.readlen).read()) * float(config["mismatch_probability"]))
        * SIGMA_MIN_SCORE_CONSTANT
    )


rule align_taxon_single_end:
    input:
        fastq="{query}/fastq/SE/{sample}_mapq.fastq.gz",
        # bt2idx="database/{orgname}/{accession}.1.bt2l",
        db_idx="database/idx_database_{query}.done",
        readlen="{query}/fastq/SE/{sample}_mapq.readlen",
    log:
        "{query}/sigma/{sample}/SE/{orgname}/{accession}.log",
    output:
        bam_file="{query}/sigma/{sample}/SE/{orgname}/{orgname}_{accession}.bam",
        bai_file="{query}/sigma/{sample}/SE/{orgname}/{orgname}_{accession}.bam.bai",
    benchmark:
        repeat(
            "benchmarks/align_taxon_single_end_{query}_{sample}_{orgname}_{accession}.benchmark.txt",
            1,
        )
    params:
        min_score=get_min_score,
        min_frag_length=MIN_FRAG_LEN,
        max_frag_length=MAX_FRAG_LEN,
    threads: config["bowtie2_treads"] # usually single threaded - the user can change it
    shell:
        "( bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
        "--score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {threads} "
        '-x database/"{wildcards.orgname}"/{wildcards.accession} '
        "-I {params.min_frag_length} -X {params.max_frag_length} "
        "-U {input.fastq} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log} "
        "; samtools index {output.bam_file}"


rule align_taxon_paired_end:
    input:
        fastq_r1="{query}/fastq/PE/{sample}_R1_mapq.fastq.gz",
        fastq_r2="{query}/fastq/PE/{sample}_R2_mapq.fastq.gz",
        # bt2idx="database/{orgname}/{accession}.1.bt2l",
        db_idx="database/idx_database_{query}.done",
        readlen="{query}/fastq/PE/{sample}_mapq_pair.readlen",
    log:
        "{query}/sigma/{sample}/PE/{orgname}/{accession}.log",
    output:
        bam_file="{query}/sigma/{sample}/PE/{orgname}/{orgname}_{accession}.bam",
        bai_file="{query}/sigma/{sample}/PE/{orgname}/{orgname}_{accession}.bam.bai",
    params:
        min_score=get_min_score,
        min_frag_length=MIN_FRAG_LEN,
        max_frag_length=MAX_FRAG_LEN,
    threads: config["bowtie2_treads"]
    shell:
        "( bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
        "--score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {threads} "
        '-x database/"{wildcards.orgname}"/{wildcards.accession} '
        "-I {params.min_frag_length} -X {params.max_frag_length} "
        "-1 {input.fastq_r1} -2 {input.fastq_r2} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log} "
        "; samtools index {output.bam_file}"

# TODO delete all of this...
# noinspection PyUnresolvedReferences
def get_bamfile_paths(wildcards):
    """
    Get all the individual bam file paths for the taxa in our database.
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

    if config["SE"] or config["PE_ANCIENT"]:
        for key, seq in sequences.iterrows():
            orgname, accession = (
                seq["species"].replace(" ", "_").replace("[", "").replace("]", ""),
                seq["GBSeq_accession-version"],
            )
            inputs.append(
                "{query}/sigma/{sample}/SE/{orgname}/{orgname}_{accession}.bam".format(
                    query=wildcards.query,
                    sample=wildcards.sample,
                    orgname=orgname,
                    accession=accession,
                )
            )

    elif config["PE_MODERN"]:
        for key, seq in sequences.iterrows():
            orgname, accession = (
                seq["species"].replace(" ", "_").replace("[", "").replace("]", ""),
                seq["GBSeq_accession-version"],
            )
            inputs.append(
                "{query}/sigma/{sample}/PE/{orgname}/{orgname}_{accession}.bam".format(
                    query=wildcards.query,
                    sample=wildcards.sample,
                    orgname=orgname,
                    accession=accession,
                )
            )

    return inputs

# TODO delete all of this...
rule all_alignments:
    input:
        get_bamfile_paths,
    output:
        "{query}/sigma/{sample}_alignments.done",
    benchmark:
        repeat("benchmarks/all_alignments_{query}_{sample}.benchmark.txt", 1)
    shell:
        "touch {output}"
