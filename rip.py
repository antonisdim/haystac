#! /usr/bin/env python
"""
Execution script for snakemake workflows.
"""
import argparse
import os.path
import snakemake
import sys
import pprint
import json
import yaml

thisdir = os.path.abspath(os.path.dirname(__file__))


def main(args):
    # first, find the Snakefile
    snakefile = os.path.join(thisdir, 'Snakefile')
    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n'.format(snakefile))
        sys.exit(-1)

    # next, find the workflow config file
    config_yaml = None
    if os.path.exists(args.config_yaml) and not os.path.isdir(args.config_yaml):
        config_yaml = args.config_yaml
    else:
        for suffix in ('', '.yaml'):
            tryfile = os.path.join(thisdir, args.workflowfile + suffix)
            if os.path.exists(tryfile) and not os.path.isdir(tryfile):
                sys.stderr.write('Found config.yaml at {}\n'.format(tryfile))
                config_yaml = tryfile
                break

    if not config_yaml:
        sys.stderr.write('Error: cannot find config.yaml {}\n'.format(args.config_yaml))
        sys.exit(-1)

    with open('config.yaml') as fin:
        config = yaml.safe_load(fin)

    print('--------')
    print('details!')
    print('\tsnakefile: {}'.format(snakefile))
    print('\tconfig: {}'.format(config_yaml))
    print('--------')

    # run!!
    status = snakemake.snakemake(snakefile, config=config, printshellcmds=True, dryrun=args.dry_run)

    if status:  # translate "success" into shell exit code of 0
        return 0
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run snakemake workflows', usage='''run <config.yaml>
    Run snakemake workflows, using the given workflow name & parameters file.
    ''')

    parser.add_argument('config_yaml')
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('--mismatch_probability', help='base mismatch probability <float> (default: 0.05)')
    parser.add_argument('--bowtie2_treads', help='threads for the bowtie2 alignments <int> (default: 1)')
    parser.add_argument('--SRA_LOOKUP', help='fetch raw data files from the SRA <bool> (default: False)')
    parser.add_argument('--PE_ANCIENT', help='treat the data as paired end end but use only COLLAPSED reads. \n'
                                             'mandatory and mutually exclusive with PE_MODERN and SE')
    parser.add_argument('--PE_MODERN', help='treat the data as paired end end end but use only COMPLETE PAIRS of '
                                            'reads. \n '
                                            'mandatory and mutually exclusive with PE_ANCIENT and SE')
    parser.add_argument('--SE', help='treat the data as single end end. \n'
                                     'mandatory and mutually exclusive with PE_ANCIENT and PE_MODERN')
    parser.add_argument('--WITH_REFSEQ_REP', help='use the prokaryotic representative species of the RefSeq DB '
                                                  'for the species id pipeline. only species no strains. '
                                                  'either or both of WITH_REFSEQ_REP and '
                                                  'WITH_ENTREZ_QUERY should be set (default: True)')
    parser.add_argument('--WITH_ENTREZ_QUERY', help='use the recordset returned from a specific entrez query in '
                                                    'the species id pipeline. written in the NCBI query language '
                                                    '(copy it from the website).either or both of WITH_REFSEQ_REP and '
                                                    'WITH_ENTREZ_QUERY should be set (default: True)')
    parser.add_argument('--SPECIFIC_GENUS', help='list containing the names of specific genera the abundances should '
                                                 'be calculated on <["genus"]>')
    parser.add_argument('--MEM_RESOURCES_MB', help='max mem resources allowed to be used ofr indexing the input for '
                                                   'the filtering alignment '
                                                   '(default: max available memory on the machine)')
    parser.add_argument('--MEM_RESCALING_FACTOR', help='factor to rescale/chunk the input file for the mutlifasta '
                                                       'index for the filtering alignemnt (default: 2.5)')

    args = parser.parse_args()

    sys.exit(main(args))
