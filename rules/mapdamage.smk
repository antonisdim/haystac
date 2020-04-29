#!/usr/bin/env python
# -*- coding: utf-8 -*-


SPECIFIC_GENUS = config['SPECIFIC_GENUS']
PE_ANCIENT = config['PE_ANCIENT']
PE_MODERN = config['PE_MODERN']
SE = config['SE']


##### Target rules #####


rule run_mapdamage:
    input:
        bam="{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam",
        ref_genome="database/{orgname}/{accession}.fasta.gz"
    log:
        "{query}/mapdamage/{sample}/{orgname}_{accession}.log"
    output:
        directory("{query}/mapdamage/{sample}/{orgname}_{accession}")
    shell:
        "mapDamage -i {input.bam} -r {input.ref_genome} -d {output}"



# noinspection PyUnresolvedReferences
def get_mapdamage_out_dir_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """
    pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
    sequences = pd.read_csv(pick_sequences.output[0], sep='\t')

    if len(sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    if WITH_REFSEQ_REP:
        refseq_rep_prok = checkpoints.entrez_refseq_accessions.get(query=wildcards.query)

        refseq_genomes = pd.read_csv(refseq_rep_prok.output[0], sep='\t')
        genbank_genomes = pd.read_csv(refseq_rep_prok.output[1], sep='\t')
        assemblies = pd.read_csv(refseq_rep_prok.output[2], sep='\t')
        refseq_plasmids = pd.read_csv(refseq_rep_prok.output[3], sep='\t')
        genbank_plasmids = pd.read_csv(refseq_rep_prok.output[4], sep='\t')

        sequences = pd.concat([sequences, refseq_genomes, genbank_genomes, assemblies, refseq_plasmids,
                               genbank_plasmids])

    inputs = []

    if SPECIFIC_GENUS:
        sequences = sequences[sequences['species'].str.contains( "|".join(SPECIFIC_GENUS))]

    for key, seq in sequences.iterrows():
        orgname, accession = seq['species'].replace(" ", "_"), seq['GBSeq_accession-version']

        inputs.append('{query}/mapdamage/{sample}/{orgname}_{accession}'.
        format(query=wildcards.query, sample=wildcards.sample, orgname=orgname, accession=accession))

    if PE_ANCIENT or SE:
        return inputs
    else:
        print("WARNING: mapDamage has not been optimised to analyse paired end alignment data.")
        return inputs
        # raise RuntimeError('PE data, mapDamage cannot run with that input format. Either collapse the reads, '
        #                    'use SE data or do not include that rule.')



rule all_mapdamage:
    input:
        get_mapdamage_out_dir_paths
    output:
        "{query}/mapdamage/{sample}_mapdamage.done"
    benchmark:
        repeat("benchmarks/all_alignments_{query}_{sample}.benchmark.txt", 3)
    shell:
        "touch {output}"