#!/usr/bin/env python
# -*- coding: utf-8 -*-

WITH_REFSEQ_REP = config['WITH_REFSEQ_REP']
WITH_ENTREZ_QUERY = config['WITH_ENTREZ_QUERY']
SRA_LOOKUP = config['SRA_LOOKUP']
PE_ANCIENT = config['PE_ANCIENT']
PE_MODERN = config['PE_MODERN']
SE = config['SE']
SPECIFIC_GENUS = config['SPECIFIC_GENUS']


##### Target rules #####


def get_inputs_for_count_fastq_len(wildcards):
    if SRA_LOOKUP:
        if PE_MODERN:
            return "fastq_inputs/PE_mod/{sample}_R1_adRm.fastq.gz".format(sample=wildcards.sample)
        elif PE_ANCIENT:
            return "fastq_inputs/PE_anc/{sample}_adRm.fastq.gz".format(sample=wildcards.sample)
        elif SE:
            return "fastq_inputs/SE/{sample}_adRm.fastq.gz".format(sample=wildcards.sample)

    else:
        if PE_MODERN:
            return config['samples'][wildcards.sample]['R1']
        elif PE_ANCIENT:
            return config['samples'][wildcards.sample]
        elif SE:
            return config['samples'][wildcards.sample]


rule count_fastq_length:
    input:
        fastq=get_inputs_for_count_fastq_len
    log:
        "fastq_inputs/meta/{sample}.log"
    output:
        "fastq_inputs/meta/{sample}.size"
    benchmark:
        repeat("benchmarks/count_fastq_length_{sample}.benchmark.txt", 3)
    shell:
        "seqtk seq -A {input.fastq} | grep -v '^>' | wc -l 1> {output} 2> {log}"



def get_bams_for_ts_tv_count(wildcards):
    if config['PE_MODERN']:
        return "{query}/sigma/{sample}/PE/{orgname}/{orgname}_{accession}.bam".format(query=wildcards.query,
            sample=wildcards.sample, orgname=wildcards.orgname, accession=wildcards.accession)
    elif config['PE_ANCIENT'] or config['SE']:
        return "{query}/sigma/{sample}/SE/{orgname}/{orgname}_{accession}.bam".format(query=wildcards.query,
            sample=wildcards.sample, orgname=wildcards.orgname, accession=wildcards.accession)



rule count_accession_ts_tv:
    input:
        get_bams_for_ts_tv_count
    output:
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}.csv"
    log:
        "{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}.log"
    benchmark:
        repeat("benchmarks/count_accession_ts_tv_{query}_{sample}_{orgname}_{accession}.benchmark.txt", 1)
    params:
        pairs=PE_MODERN
    script:
        "../scripts/count_accession_ts_tv.py"


# noinspection PyUnresolvedReferences
def get_ts_tv_count_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """

    sequences = pd.DataFrame()

    if WITH_ENTREZ_QUERY:
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

        if WITH_ENTREZ_QUERY:
            sequences = pd.concat([sequences, refseq_genomes, genbank_genomes, assemblies, refseq_plasmids,
                                   genbank_plasmids])
        else:
            sequences = pd.concat([refseq_genomes, genbank_genomes, assemblies, refseq_plasmids,
                                   genbank_plasmids])

    inputs = []

    if SPECIFIC_GENUS:
        sequences = sequences[sequences['species'].str.contains( "|".join(SPECIFIC_GENUS))]

    for key, seq in sequences.iterrows():
        orgname, accession = seq['species'].replace(" ", "_"), seq['GBSeq_accession-version']

        inputs.append('{query}/ts_tv_counts/{sample}/{orgname}_count_{accession}.csv'.
        format(query=wildcards.query, sample=wildcards.sample, orgname=orgname, accession=accession))

    return inputs


rule initial_ts_tv:
    input:
        get_ts_tv_count_paths
    output:
        "{query}/ts_tv_counts/{sample}/all_ts_tv_counts.csv"
    log:
        "{query}/ts_tv_counts/{sample}/all_ts_tv_counts.log"
    benchmark:
        repeat("benchmarks/initial_ts_tv_{query}_{sample}.benchmark.txt", 1)
    shell:
        "cat {input} 1> {output} 2> {log}"


def get_right_readlen(wildcards):
    if PE_MODERN:
        return "{query}/fastq/PE/{sample}_mapq_pair.readlen".format(query=wildcards.query, sample=wildcards.sample)
    else:
        return "{query}/fastq/SE/{sample}_mapq.readlen".format(query=wildcards.query, sample=wildcards.sample)


rule calculate_likelihoods:
    input:
        "{query}/ts_tv_counts/{sample}/all_ts_tv_counts.csv",
        get_right_readlen,
        get_ts_tv_count_paths
    output:
        "{query}/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
        "{query}/probabilities/{sample}/{sample}_probability_model_params.json"
    log:
        "{query}/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.log"
    benchmark:
        repeat("benchmarks/calculate_likelihoods_{query}_{sample}.benchmark.txt", 1)
    script:
        "../scripts/calculate_likelihoods.py"  #todo ask Evan to check if they are the same with the SQL commands



rule calculate_taxa_probabilities:
    input:
        "{query}/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
        "{query}/probabilities/{sample}/{sample}_probability_model_params.json",
        "fastq_inputs/meta/{sample}.size"
    output:
        "{query}/probabilities/{sample}/{sample}_posterior_probabilities.csv"
    log:
        "{query}/probabilities/{sample}/{sample}_posterior_probabilities.log"
    benchmark:
        repeat("benchmarks/calculate_taxa_probabilities_{query}_{sample}.benchmark.txt", 1)
    params:
        submatrices=False
    script:
        "../scripts/calculate_taxa_probabilities.py"



rule fasta_idx:
    input:
        "database/{orgname}/{accession}.fasta.gz"
    output:
        "database/{orgname}/{accession}.fasta.gz.fai"
    log:
        "database/{orgname}/{accession}.fasta.gz.fai.log"
    benchmark:
        repeat("benchmarks/fasta_idx_{orgname}_{accession}.benchmark.txt", 3)
    shell:
        "samtools faidx {input} 2> {log}"



rule coverage_t_test:
    input:
        "{query}/sigma/{sample}/{reads}/{orgname}/{orgname}_{accession}.bam",
        "database/{orgname}/{accession}.fasta.gz.fai"
    output:
        "{query}/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}_{reads}.txt"
    log:
        "{query}/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}_{reads}.log"
    benchmark:
        repeat("benchmarks/coverage_t_test_{query}_{sample}_{orgname}_{accession}_{reads}.benchmark.txt", 1)
    script:
        "../scripts/coverage_t_test.py"



# noinspection PyUnresolvedReferences
def get_t_test_values_paths(wildcards):
    """
    Get all the individual cav file paths for the taxa in our database.
    """

    sequences = pd.DataFrame()

    if WITH_ENTREZ_QUERY:
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

        if WITH_ENTREZ_QUERY:
            sequences = pd.concat([sequences, refseq_genomes, genbank_genomes, assemblies, refseq_plasmids,
                                   genbank_plasmids])
        else:
            sequences = pd.concat([refseq_genomes, genbank_genomes, assemblies, refseq_plasmids,
                                   genbank_plasmids])

    inputs = []

    if SPECIFIC_GENUS:
        sequences = sequences[sequences['species'].str.contains( "|".join(SPECIFIC_GENUS))]

    reads = ''
    if PE_MODERN:
        reads = 'PE'
    elif PE_ANCIENT or SE:
        reads = 'SE'

    for key, seq in sequences.iterrows():
        orgname, accession = seq['species'].replace(" ", "_"), seq['GBSeq_accession-version']

        inputs.append('{query}/probabilities/{sample}/{orgname}_t_test_pvalue_{accession}_{reads}.txt'.
        format(query=wildcards.query, sample=wildcards.sample, orgname=orgname, accession=accession, reads=reads))

    return inputs


rule cat_pvalues:
    input:
        get_t_test_values_paths
    output:
        "{query}/probabilities/{sample}/{sample}_t_test_pvalues.txt",
    log:
        "{query}/probabilities/{sample}/{sample}_t_test_pvalues.log"
    benchmark:
        repeat("benchmarks/cat_pvalues_{query}_{sample}.benchmark.txt", 3)
    shell:
        "cat {input} 1> {output} 2> {log}"



rule calculate_dirichlet_abundances:
    input:
        "{query}/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv",
        "{query}/probabilities/{sample}/{sample}_t_test_pvalues.txt",
        "fastq_inputs/meta/{sample}.size"
    output:
        "{query}/probabilities/{sample}/{sample}_posterior_abundance.tsv"
    log:
        "{query}/probabilities/{sample}/{sample}_posterior_abundance.log"
    benchmark:
        repeat("benchmarks/calculate_dirichlet_abundances_{query}_{sample}.benchmark.txt", 1)
    script:
        "../scripts/calculate_dirichlet_abundances.py"
