#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd


def entrez_refseq_create_files(input_file, nuccore_genomes_out, genbank_genomes_out, assemblies_out,
                               nuccore_plasmids_out, genbank_plasmids_out):
    prok_refseq_rep = pd.read_csv(input_file, sep='\t')

    prok_refseq_rep_rmdup = prok_refseq_rep[~prok_refseq_rep['#Species/genus'].duplicated()]

    assemblies = prok_refseq_rep_rmdup.loc[prok_refseq_rep_rmdup['WGS'].notna(), ['#Species/genus', 'WGS']]

    assemblies['#Species/genus'] = assemblies['#Species/genus'].str.replace(' ', '_')

    nuccore = prok_refseq_rep_rmdup.loc[
        prok_refseq_rep_rmdup['Chromosome RefSeq'].notna(), ['#Species/genus', 'Chromosome RefSeq']]

    nuccore['#Species/genus'] = nuccore['#Species/genus'].str.replace(' ', '_')

    genbank = prok_refseq_rep_rmdup.loc[
        prok_refseq_rep_rmdup['Chromosome GenBank'].notna(), ['#Species/genus', 'Chromosome GenBank']]

    genbank['#Species/genus'] = genbank['#Species/genus'].str.replace(' ', '_')

    genbank_filtered = genbank[(~genbank['#Species/genus'].isin(assemblies['#Species/genus'])) & (
        ~genbank['#Species/genus'].isin(nuccore['#Species/genus']))]

    nuccore_filtered = nuccore[(~nuccore['#Species/genus'].isin(assemblies['#Species/genus'])) & (
        ~nuccore['#Species/genus'].isin(genbank_filtered['#Species/genus']))]

    assemblies_filtered = assemblies[(~assemblies['#Species/genus'].isin(nuccore_filtered['#Species/genus'])) & (
        ~assemblies['#Species/genus'].isin(genbank_filtered['#Species/genus']))]

    nuccore_plasmids = prok_refseq_rep_rmdup[
                           prok_refseq_rep_rmdup['Plasmid RefSeq'].notna() &
                           prok_refseq_rep_rmdup['WGS'].isna()].loc[:, ['#Species/genus', 'Plasmid RefSeq']]

    genbank_plasmids = prok_refseq_rep_rmdup[
                           prok_refseq_rep_rmdup['Plasmid GenBank'].notna() & prok_refseq_rep_rmdup['WGS'].isna()].loc[
                       :, ['#Species/genus', 'Plasmid GenBank']]
    genbank_plasmids_filtered = genbank_plasmids[
        ~genbank_plasmids['#Species/genus'].isin(nuccore_plasmids['#Species/genus'])]

    nuccore_plasmids.loc[:, 'Plasmid RefSeq'] = nuccore_plasmids['Plasmid RefSeq'].str.split(',')
    nuccore_plasmids_exploded = nuccore_plasmids.explode('Plasmid RefSeq')

    # todo these lines give me this warning: SettingWithCopyWarning:
    #  A value is trying to be set on a copy of a slice from a DataFrame.
    #  Read about it, but can't figure out why it's happening
    genbank_plasmids_filtered.loc[:, 'Plasmid GenBank'] = genbank_plasmids_filtered['Plasmid GenBank'].str.split(',')
    genbank_plasmids_filtered_exploded = genbank_plasmids_filtered.explode('Plasmid GenBank')

    header = ['species', 'GBSeq_accession-version']

    genbank_filtered.to_csv(genbank_genomes_out, sep='\t', header=header, index=False)
    nuccore_filtered.to_csv(nuccore_genomes_out, sep='\t', header=header, index=False)
    assemblies_filtered.to_csv(assemblies_out, sep='\t', header=header, index=False)
    genbank_plasmids_filtered_exploded.to_csv(genbank_plasmids_out, sep='\t', header=header, index=False)
    nuccore_plasmids_exploded.to_csv(nuccore_plasmids_out, sep='\t', header=header, index=False)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_refseq_create_files(
        input_file=snakemake.input[0],
        nuccore_genomes_out=snakemake.output[0],
        genbank_genomes_out=snakemake.output[1],
        assemblies_out=snakemake.output[2],
        nuccore_plasmids_out=snakemake.output[3],
        genbank_plasmids_out=snakemake.output[4]
    )
