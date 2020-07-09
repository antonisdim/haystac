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
    snakefile = os.path.join(thisdir, "Snakefile")
    if not os.path.exists(snakefile):
        sys.stderr.write("Error: cannot find Snakefile at {}\n".format(snakefile))
        sys.exit(-1)

    # next, find the workflow config file

    if os.path.exists(args.config_yaml):
        with open(args.config_yaml) as fin:
            config_yaml = yaml.safe_load(fin)
    else:
        with open(os.path.join(thisdir, "config.yaml")) as fin:
            config_yaml = yaml.safe_load(fin)

    if not config_yaml:
        sys.stderr.write(
            "Error: cannot find the config yaml file {}. You definitely need one, even if it's empty "
            "if you want to run everything from the CLI flags\n".format(
                args.config_yaml
            )
        )
        sys.exit(-1)

    config_args = vars(args)

    # config = {**config_yaml, **config_args}
    config_yaml.update((k, v) for k, v in config_args.items() if v is not None)
    config = config_yaml

    input_mode = ""
    if config["PE_ANCIENT"]:
        input_mode = "PE_anc"
    elif config["PE_MODERN"]:
        input_mode = "PE_mod"
    elif config["SE"]:
        input_mode = "SE"

    data_preprocessing = ""
    if config["PE_MODERN"]:
        data_preprocessing = "fastq_inputs/{input_mode}/{sample}_R1_adRm.fastq.gz".format(
            sample=config["sample_name"], input_mode=input_mode
        )
    elif config["PE_ANCIENT"] or config["SE"]:
        data_preprocessing = "fastq_inputs/{input_mode}/{sample}_adRm.fastq.gz".format(
            sample=config["sample_name"], input_mode=input_mode
        )

    entrez_build_prok_refseq_rep = "{query}/bowtie/refseq_rep_refseq_prok.fasta.gz".format(
        query=config["query_name"]
    )
    entrez = "{query}/bowtie/{query}_entrez.fasta.gz".format(query=config["query_name"])

    bowtie = ""
    if config["PE_MODERN"]:
        bowtie = "{query}/fastq/PE/{sample}_mapq_pair.readlen".format(
            query=config["query_name"], sample=config["sample_name"]
        )
    elif config["PE_ANCIENT"] or config["SE"]:
        bowtie = "{query}/fastq/SE/{sample}_mapq.readlen".format(
            query=config["query_name"], sample=config["sample_name"]
        )

    index_database = "database/idx_database_{query}.done".format(
        query=config["query_name"]
    )

    bt2_filter_idx = "{query}/bowtie/bowtie_index.done".format(
        query=config["query_name"]
    )
    bowtie_meta = "{query}/sigma/{sample}_alignments.done".format(
        query=config["query_name"], sample=config["sample_name"]
    )

    metagenomics_probabilities = "{query}/probabilities/{sample}/{sample}_posterior_probabilities.csv".format(
        query=config["query_name"], sample=config["sample_name"]
    )
    metagenomics_abundances = "{query}/probabilities/{sample}/{sample}_posterior_abundance.tsv".format(
        query=config["query_name"], sample=config["sample_name"]
    )

    mapdamage = "{query}/mapdamage/{sample}_mapdamage.done".format(
        query=config["query_name"], sample=config["sample_name"]
    )

    target_list = [
        bowtie,
        bowtie_meta,
        metagenomics_probabilities,
        metagenomics_abundances,
        mapdamage,
    ]

    if config["WITH_DATA_PREPROCESSING"]:
        target_list.append(data_preprocessing)

    if config["WITH_ENTREZ_QUERY"]:
        target_list.append(entrez)
    if config["WITH_REFSEQ_REP"]:
        target_list.append(entrez_build_prok_refseq_rep)

    if args.data_preprocess:
        target_list = [data_preprocessing]

    if args.build_db:
        target_list = []
        if config["WITH_ENTREZ_QUERY"]:
            target_list.append(entrez)
        if config["WITH_REFSEQ_REP"]:
            target_list.append(entrez_build_prok_refseq_rep)

    if args.index_db:
        target_list = [index_database]

    if args.aln_filter_index:
        target_list = [bt2_filter_idx]

    if args.aln_filter:
        target_list = [bowtie]

    if args.aln_meta:
        target_list = [bowtie_meta]

    if args.probabilities:
        target_list = [metagenomics_probabilities]

    if args.abundances:
        target_list = [metagenomics_abundances]

    if args.mapdamage:
        target_list = [mapdamage]

    print("--------")
    print("details!")
    print("\tsnakefile: {}".format(snakefile))
    print("\tconfig: {}".format(config))
    print("\ttargets: {}".format(target_list))
    print("--------")

    # run!!
    status = snakemake.snakemake(
        snakefile,
        config=config,
        targets=target_list,
        printshellcmds=True,
        dryrun=args.dry_run,
        cores=int(args.cores),
        keepgoing=True,
        restart_times=15,
        touch=args.touch,
        forceall=args.touch,
        lock=args.no_lock,
        unlock=args.unlock,
    )

    if status:  # translate "success" into shell exit code of 0
        return 0
    return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="run snakemake workflows",
        usage="""run <config.yaml>
    Run snakemake workflows, using the given workflow name & parameters file.
    """,
    )

    parser.add_argument("config_yaml")
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument(
        "--query_name",
        help="name of the query and the output directory <str>",
        metavar="",
    )
    parser.add_argument(
        "--entrez_query",
        help="Actual NCBI query in the NCBI query language",
        metavar="",
    )
    parser.add_argument(
        "--entrez_email", help="email address for NCBI identification", metavar=""
    )
    parser.add_argument(
        "--entrez_batchsize",
        help="batchsize for fetching records from NCBI",
        metavar="",
    )
    parser.add_argument(
        "--entrez_rank",
        help="Taxonomic rank to perform the identifications on",
        metavar="",
    )
    parser.add_argument("--sample_name", help="Sample name prefix", metavar="")
    parser.add_argument(
        "--sample_fastq",
        help="Path to the fastq input file. Can be raw or with adapters " "removed",
        metavar="",
    )
    parser.add_argument(
        "--sample_fastq_R1",
        help="Path to the mate 1 fastq input file, if reads are PE. "
        "Can be raw or with adapters removed",
        metavar="",
    )
    parser.add_argument(
        "--sample_fastq_R2",
        help="Path to the mate 2 fastq input file, if reads are PE. "
        "Can be raw or with adapters removed",
        metavar="",
    )
    parser.add_argument(
        "--mismatch_probability",
        help="base mismatch probability <float> (default: 0.05)",
        metavar="",
    )
    parser.add_argument(
        "--bowtie2_treads",
        help="threads for the bowtie2 alignments <int> (default: 1)",
        metavar="",
    )
    parser.add_argument(
        "--SRA_LOOKUP",
        help="fetch raw data files from the SRA <bool> (default: False)",
        metavar="",
    )
    parser.add_argument(
        "--PE_ANCIENT",
        help="treat the data as paired end end but use only COLLAPSED reads. \n"
        "mandatory and mutually exclusive with PE_MODERN and SE",
        metavar="",
    )
    parser.add_argument(
        "--PE_MODERN",
        help="treat the data as paired end end end but use only COMPLETE PAIRS of "
        "reads. \n "
        "mandatory and mutually exclusive with PE_ANCIENT and SE",
        metavar="",
    )
    parser.add_argument(
        "--SE",
        help="treat the data as single end end. \n"
        "mandatory and mutually exclusive with PE_ANCIENT and PE_MODERN",
        metavar="",
    )
    parser.add_argument(
        "--WITH_REFSEQ_REP",
        help="use the prokaryotic representative species of the RefSeq DB "
        "for the species id pipeline. only species no strains. "
        "either or both of WITH_REFSEQ_REP and "
        "WITH_ENTREZ_QUERY should be set (default: True)",
        metavar="",
    )
    parser.add_argument(
        "--WITH_ENTREZ_QUERY",
        help="use the recordset returned from a specific entrez query in "
        "the species id pipeline. written in the NCBI query language "
        "(copy it from the website).either or both of WITH_REFSEQ_REP and "
        "WITH_ENTREZ_QUERY should be set (default: True)",
        metavar="",
    )
    parser.add_argument(
        "--WITH_DATA_PREPROCESSING",
        help="Remove adapters from raw fastq files",
        metavar="",
    )
    parser.add_argument(
        "--SPECIFIC_GENUS",
        nargs="+",
        help="list containing the names of specific genera "
        "the abundances should be calculated "
        'on <["genus"]>',
        metavar="",
    )
    parser.add_argument(
        "--MEM_RESOURCES_MB",
        help="max mem resources allowed to be used ofr indexing the input for "
        "the filtering alignment "
        "(default: max available memory on the machine)",
        metavar="",
    )
    parser.add_argument(
        "--MEM_RESCALING_FACTOR",
        help="factor to rescale/chunk the input file for the mutlifasta "
        "index for the filtering alignment (default: 2.5)",
        metavar="",
    )
    parser.add_argument(
        "--data_preprocess",
        action="store_true",
        help="Only process raw fastq files, "
        "incl downloading fastq files from the SRA, "
        "compressing them and removing any sequencing "
        "adapters  <bool> (default: False)",
    )
    parser.add_argument(
        "--build_db",
        action="store_true",
        help="ONly build the database from refseq and/or "
        "other sequences retrieved from an NCBI query",
    )
    parser.add_argument(
        "--index_db",
        action="store_true",
        help="Index each taxon's reference genome that is "
        "contained in the reference database, using bowtie2 "
        " <bool> (default: False)",
    )
    parser.add_argument(
        "--aln_filter_index",
        action="store_true",
        help="Build the index required for the filtering "
        "alignment step, out of all the sequences "
        "contained in the reference database "
        " <bool> (default: False)",
    )
    parser.add_argument(
        "--aln_filter",
        action="store_true",
        help="Perform the filtering alignment step " " <bool> (default: False)",
    )
    parser.add_argument(
        "--aln_meta",
        action="store_true",
        help="Perform metagenomic alignments. "
        "The aligned reads from the filtering alignment step "
        "are aligned to each reference genomes separately "
        "using bowtie2  <bool> (default: False)",
    )
    parser.add_argument(
        "--probabilities",
        action="store_true",
        help="Taxon assignment posterior probabilities are "
        "calculated for a given sample  <bool> "
        "(default: False)",
    )
    parser.add_argument(
        "--abundances",
        action="store_true",
        help="Mean posterior abundances for each taxon in our "
        "ref db, for a given sample, are calculated "
        " <bool> (default: False)",
    )
    parser.add_argument(
        "--mapdamage",
        action="store_true",
        help="Performs a chemical damage analysis, "
        "by using mapDamage  <bool> (default: False)",
    )
    parser.add_argument(
        "--touch",
        action="store_true",
        help="Touches the timestamps of the requested target files "
        "AND ALL of their dependencies  <bool> (default: False)",
    )
    parser.add_argument(
        "-j", "--cores", help="Number of cores for smk to use", metavar=""
    )
    parser.add_argument(
        "-l",
        "--no_lock",
        action="store_true",
        help="Apply no lock to the working directory " " <bool> (default: False)",
    )
    parser.add_argument(
        "-ul",
        "--unlock",
        action="store_true",
        help="Unlock the working directory after smk is "
        "abruptly killed  <bool> (default: False)",
    )

    args = parser.parse_args()

    sys.exit(main(args))
