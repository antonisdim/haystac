#! /usr/bin/env python
"""
Execution script for snakemake workflows.
"""
__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import datetime

import pathlib

import argcomplete
import argparse
import os.path
import re
import sys
import os

import snakemake
import yaml

from multiprocessing import cpu_count
from Bio import Entrez
from psutil import virtual_memory
from pathlib import Path

from workflow.scripts.utilities import (
    ValidationError,
    EmailType,
    WritablePathType,
    PositiveIntType,
    FloatRangeType,
    IntRangeType,
    BoolType,
    JsonType,
)

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

MAX_CPU = cpu_count()
MAX_MEM_MB = int(virtual_memory().total / 1024 ** 2)

CONFIG_DEFAULT = "./config/config.yaml"
CONFIG_USER = pathlib.Path("~/.haystack/config.yaml").expanduser()

DATABASE_MODES = ["fetch", "index", "build"]
TAXONOMIC_RANKS = ["genus", "species", "subspecies", "serotype"]

RESTART_TIMES = 3

BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# maximum concurrent Entrez requests
MAX_ENTREZ_REQUESTS = 3


class Haystack(object):
    def __init__(self):
        parser = argparse.ArgumentParser(
            description="HAYSTACK: A Bayesian framework for robust and rapid species identification",
            usage="""haystack <command> [<args>]

The haystack commands are:
   config         Advanced configuration options for haystack
   database       Build a database of target species
   sample         Prepare a sample for analysis
   analyse        Analyse a sample against a database
   
""",
        )
        parser.add_argument(
            "command", choices=["config", "database", "sample", "analyse"], help="Command to run",
        )

        # get the command
        argcomplete.autocomplete(parser)  # TODO autocomplete does not work
        args = parser.parse_args(sys.argv[1:2])

        # load the default config
        with open(CONFIG_DEFAULT) as fin:
            self.config_default = yaml.safe_load(fin)

            # resolve the home directory of the user (i.e. ~/)
            self.config_default["cache"] = str(pathlib.Path(self.config_default["cache"]).expanduser())

        try:
            # load the user config
            with open(CONFIG_USER) as fin:
                self.config_user = yaml.safe_load(fin)

        except FileNotFoundError:
            # make the config folder (if necessary)
            os.makedirs(os.path.dirname(CONFIG_USER), exist_ok=True)

            self.config_user = dict()

        if args.command != "config" and not self.config_user.get("email"):
            # email address is mandatory
            print(
                "Before using haystack, please set your email address. This is a required step for running NCBI "
                "queries\n\n`haystack config --email <address>`"
            )
            exit(1)

        # merge the config dictionaries (giving precedence to the config_user)
        self.config_merged = z = {**self.config_default, **self.config_user}

        try:
            # use dispatch pattern to invoke method with same name
            getattr(self, args.command)()

        except ValidationError as error:
            print(f"haystack: error: {error}")
            exit(1)

    def config(self):
        """
        Advanced configuration options for haystack
        """
        parser = argparse.ArgumentParser(description="Advanced configuration options for haystack")

        parser.add_argument(
            "--email",
            help="Email address for NCBI identification (mandatory).",
            metavar="<address>",
            type=EmailType("RFC5322"),
        )

        parser.add_argument(
            "--cache",
            help=f"Cache folder for all the genomes downloaded from NCBI (default: {self.config_default['cache']}).",
            metavar="<path>",
            type=WritablePathType(),
            default=self.config_default["cache"],
        )

        parser.add_argument(
            "--batchsize",
            help=f"Batch size for fetching records from NCBI (default: {self.config_default['batchsize']})",
            type=PositiveIntType(),
            metavar="<int>",
            default=self.config_default["batchsize"],
        )

        parser.add_argument(
            "--mismatch-probability",
            help=f"Base mismatch probability (default: {self.config_default['mismatch_probability']})",
            type=FloatRangeType(0.01, 0.10),
            metavar="<float>",
            default=self.config_default["mismatch_probability"],
        )

        parser.add_argument(
            "--bowtie2-threads",
            help=f"Number of threads to use for each bowtie2 alignment "
            f"(default: {self.config_default['bowtie2_threads']})",
            type=IntRangeType(1, MAX_CPU),
            metavar="<int>",
            default=self.config_default["bowtie2_threads"],
        )

        parser.add_argument(
            "--bowtie2-scaling",
            help=f"Rescaling factor to keep the bowtie2 mutlifasta index below the maximum memory limit "
            f"(default: {self.config_default['bowtie2_scaling']})",
            type=FloatRangeType(0, 100),
            metavar="<float>",
            default=self.config_default["bowtie2_scaling"],
        )

        parser.add_argument(
            "--use-conda",
            help=f"Use conda as a package manger (default: {self.config_default['use_conda']})",
            type=BoolType(),
            metavar="<bool>",
            default=self.config_default["use_conda"],
        )

        # now that we're inside a subcommand, ignore the first two arguments
        argcomplete.autocomplete(parser)
        args = parser.parse_args(sys.argv[2:])

        config_user = dict()

        # get the user choices that differ from the defaults
        for key, value in vars(args).items():
            if value != self.config_default.get(key) and value is not None:
                config_user[key] = value

        # save the user config
        with open(CONFIG_USER, "w") as fout:
            yaml.safe_dump(config_user, fout, default_flow_style=False)

    def _common_arguments(self, parser):
        """
        Add the common arguments shared by the `database`, `sample` and `analyse` commands
        """
        parser.add_argument(
            "--cores",
            help=f"Maximum number of CPU cores to use (default: {self.config_default['cores']}).",
            metavar="<int>",
            type=IntRangeType(1, MAX_CPU),
            default=MAX_CPU if self.config_default["cores"] == "all" else self.config_default["cores"],
        )

        parser.add_argument(
            "--mem",
            help=f"Maximum megabytes of memory to use (default: {self.config_default['mem']}).",
            type=IntRangeType(1024, MAX_MEM_MB),
            default=MAX_MEM_MB if self.config_default["mem"] == "all" else self.config_default["mem"],
            metavar="<int>",
        )

        parser.add_argument(
            "--unlock",
            help="Unlock the output directory following a crash or hard restart (default: False).",
            action="store_true",
        )

        parser.add_argument(
            "--debug", help="Enable debugging mode (default: False)", action="store_true",
        )

        parser.add_argument(
            "--snakemake",
            help="Pass additional flags to the `snakemake` scheduler.",
            metavar="'<json>'",
            type=JsonType(),
        )

    def database(self):
        """
        Build a database of target species
        """
        parser = argparse.ArgumentParser(description="Build a database of target species")

        parser.add_argument(
            "--mode",
            choices=DATABASE_MODES,
            help=f"Database creation mode for haystack (default: {self.config_default['mode']}).\n"
            f"Alternatively choose 'fetch' to download the genomes, then 'index' to build the alignment indices.",
            metavar="<mode>",
            default=self.config_default["mode"],
        )

        parser.add_argument(
            "--rank",
            help=f"Taxonomic rank to perform the identifications on [{', '.join(TAXONOMIC_RANKS)}] "
            f"(default: {self.config_default['rank']})",
            choices=TAXONOMIC_RANKS,
            default=self.config_default["rank"],
            metavar="<rank>",
        )

        parser.add_argument(
            "--refseq-rep",
            help="Include all prokaryotic species (no strains) from the representative RefSeq DB (default: False)",
            action="store_true",
        )

        # TODO how can we validate these queries?
        parser.add_argument(
            "--query",
            help="Database query in the NCBI query language. "
            "Please refer to the documentation for assistance with constructing a valid query.",
            metavar="<query>",
        )

        parser.add_argument(
            "--query-file",
            help="File containing a database query in the NCBI query language.",
            metavar="<path>",
            type=argparse.FileType("r"),
        )

        parser.add_argument(
            "--mtDNA",
            help="For eukaryotes, download mitochondrial genomes only. "
            "Not to be used with --refseq-rep or queries containing prokaryotes (default: False)",
            action="store_true",
        )

        # TODO validate that this is 2-column and tab delimited
        parser.add_argument(
            "--accessions",
            help="Tab delimited file containing one record per row: the name of the taxon, "
            "and a valid NCBI accession code from the nucleotide, assembly or WGS databases.",
            metavar="<path>",
            type=argparse.FileType("r"),
        )

        # TODO validate that this is 3-column and tab delimited
        parser.add_argument(
            "--sequences",
            help="Tab delimited file containing one record per row: the name of the taxon, a user defined "
            "accession code, and the path to the fasta file (optionally compressed).",
            metavar="<path>",
            type=argparse.FileType("r"),
        )

        parser.add_argument(
            "--seed",
            help=f"Random seed for database indexing (default: {self.config_default['seed']})",
            metavar="<int>",
            type=int,
            default=self.config_default["seed"],
        )

        parser.add_argument(
            "--genera",
            help="Optional list of genera to restrict the abundance calculations.",
            metavar="<genus>",
            nargs="+",
            default=[],
        )

        # add the common arguments
        self._common_arguments(parser)

        parser.add_argument(
            "--output",
            help="Path to the database output directory (mandatory).",
            metavar="<path>",
            dest="db_output",
            type=WritablePathType(),
            required=True,
        )

        # print the help
        if len(sys.argv) == 2:
            parser.print_help()
            parser.exit()

        # now that we're inside a subcommand, ignore the first two arguments
        argcomplete.autocomplete(parser)
        args = parser.parse_args(sys.argv[2:])

        # must specify at least one source for the database
        if not (args.refseq_rep or args.query or args.query_file or args.accessions or args.sequences):
            raise ValidationError(
                "Please specify at least one of --refseq-rep, --query, --query-file, --accessions or --sequences"
            )

        if args.mtDNA and args.refseq_rep:
            raise ValidationError("Please specify either `--mtDNA` or `--refseq-rep` but not both.")

        if args.query and args.query_file:
            raise ValidationError("Please specify either `--query <query>` or `--query-file <path>` but not both.")

        if args.query_file:
            # load the query file
            args.query = args.query_file.read().strip()

            if not args.query:
                raise ValidationError(f"The query file '{args.query_file.name}' is empty.")

        # resolve relative paths
        args.db_output = os.path.abspath(args.db_output)

        # add all command line options to the merged config
        config = {**self.config_merged, **vars(args)}

        config_fetch = os.path.join(args.db_output, "database_fetch_config.yaml")
        config_build = os.path.join(args.db_output, "database_build_config.yaml")

        target_list = []

        if args.mode == "fetch":
            # TODO can this be replaced with a single target? and the conditional logic moved into an input function
            if args.query:
                target_list.append("bowtie/entrez_query.fasta.gz")
            if args.refseq_rep:
                target_list.append("bowtie/refseq_prok.fasta.gz")
            if args.sequences:
                target_list.append("bowtie/custom_seqs.fasta.gz")
            if args.accessions:
                target_list.append("bowtie/custom_acc.fasta.gz")

            with open(config_fetch, "w") as fout:
                yaml.safe_dump(config, fout, default_flow_style=False)

            print("Please run `haystack database --mode index` after this step.")

        elif args.mode == "index":
            target_list.append("bowtie/bowtie_index.done")

            try:
                with open(config_fetch, "r") as fin:
                    config = yaml.safe_load(fin)
            except FileNotFoundError:
                raise ValidationError(
                    "Please run haystack `database --mode fetch` before attempting to index the database."
                )

        elif args.mode == "build":
            target_list.append("idx_database.done")
            target_list.append("bowtie/bowtie_index.done")

            if os.path.exists(config_fetch):
                raise ValidationError(
                    "Please run haystack `database --mode index` as the database has alrady been fetched."
                )

            with open(config_build, "w") as fout:
                yaml.safe_dump(config, fout, default_flow_style=False)

        config["workflow_dir"] = os.path.join(BASE_DIR, "workflow")  # TODO tidy up
        config["mtDNA"] = str(args.mtDNA).lower()
        config["sequences"] = config["sequences"] or ""

        target_list = [os.path.join(args.db_output, target) for target in target_list]
        snakefile = os.path.join(BASE_DIR, "workflow/database.smk")

        return self._run_snakemake(snakefile, args, config, target_list)

    def sample(self):
        """
        Prepare a sample for analysis
        """
        parser = argparse.ArgumentParser(description="Prepare a sample for analysis")

        # TODO do we need this?
        parser.add_argument(
            "--sample-prefix", help="Sample prefix for all the future analysis.", metavar="<prefix>",
        )

        parser.add_argument(
            "--fastq",
            help="Single-end fastq input file (optionally compressed).",
            metavar="<path>",
            type=argparse.FileType("r"),
        )

        parser.add_argument(
            "--fastq-r1",
            help="Paired-end forward strand (R1) fastq input file.",
            metavar="<path>",
            type=argparse.FileType("r"),
        )

        parser.add_argument(
            "--fastq-r2",
            help="Paired-end reverse strand (R2) fastq input file.",
            metavar="<path>",
            type=argparse.FileType("r"),
        )

        parser.add_argument(
            "--collapse",
            help=f"Collapse overlapping paired-end reads, e.g. for aDNA (default: {self.config_default['collapse']})",
            type=BoolType(),
            metavar="<bool>",
            default=self.config_default["collapse"],
        )

        # TODO how to validate these accessions? do they follow a standard pattern?
        parser.add_argument(
            "--sra",
            help="Download fastq input from the SRA database",
            metavar="<accession>",
        )

        parser.add_argument(
            "--trim-adapters",
            help=f"Automatically trim sequencing adapters from fastq input "
            f"(default: {self.config_default['trim_adapters']})",
            type=BoolType(),
            metavar="<bool>",
            default=self.config_default["trim_adapters"],
        )

        # TODO does this work?
        parser.add_argument(
            "--adaperremoval-flags",
            help="Pass additional flags to `AdapterRemoval`.",
            metavar="'<json>'",
            type=JsonType(),
        )

        # add the common arguments
        self._common_arguments(parser)

        parser.add_argument(
            "--output",
            help="Path to the sample output directory (mandatory).",
            metavar="<path>",
            dest="sample_output_dir",
            type=WritablePathType(),
            required=True,
        )

        # print the help
        if len(sys.argv) == 2:
            parser.print_help()
            parser.exit()

        # now that we're inside a subcommand, ignore the first two arguments
        argcomplete.autocomplete(parser)
        args = parser.parse_args(sys.argv[2:])

        # must specify exactly one source for the sample
        if bool(args.fastq) + (bool(args.fastq_r1) and bool(args.fastq_r2)) + bool(args.sra) != 1:
            raise ValidationError(
                "Please specify either --sra or --fastq, or both --fastq-r1 and --fastq-r2."
            )

        if args.collapse and not (args.fastq_r1 and args.fastq_r2):
            raise ValidationError(
                "Collapse can only be used with --fastq-r1 and --fastq-r2."
            )

        if not args.sample_prefix and not args.sra:
            raise ValidationError(
                "Please provide a --sample-prefix name."
            )

        # resolve relative paths
        args.sample_output_dir = os.path.abspath(args.sample_output_dir)

        # cast all paths as strings
        args.fastq = args.fastq.name if args.fastq else ""
        args.fastq_r1 = args.fastq_r1.name if args.fastq_r1 else ""
        args.fastq_r2 = args.fastq_r2.name if args.fastq_r2 else ""

        # add all command line options to the merged config
        config = {**self.config_merged, **vars(args)}

        # TODO are these flags really necessary?
        #      why not combine into one single config value? (e.g. config["mode"] = SE|PE|COLLAPSE
        config["SE"] = not (args.fastq_r1 and args.fastq_r2)
        config["PE_ANCIENT"] = not config["SE"] and args.collapse
        config["PE_MODERN"] = not config["SE"] and not args.collapse

        # TODO rethink this...
        if args.sra:
            if args.sample_prefix:
                raise ValidationError(
                    "--sample-prefix cannot be used with and SRA accession."
                )

            # use the SRA accession as the sample prefix
            config['sample_prefix'] = args.sra

            # get paired/single status of the accession
            Entrez.email = config["email"]
            sra_id = Entrez.read(Entrez.esearch(db="sra", term=config["sra"]))["IdList"]
            if "paired" in str(Entrez.read(Entrez.esummary(db="sra", id=sra_id))).lower():
                if config["collapse"]:
                    config["PE_ANCIENT"] = True
                else:
                    config["PE_MODERN"] = True

                config["fastq_r1"] = f"sra_data/PE/{config['sra']}_R1.fastq.gz"
                config["fastq_r2"] = f"sra_data/PE/{config['sra']}_R1.fastq.gz"

            else:
                config["SE"] = True
                config["fastq"] = f"sra_data/SE/{config['sra']}.fastq.gz"

        target_list = []
        target_list.append(f"fastq_inputs/meta/{config['sample_prefix']}.size")

        if config["trim_adapters"]:
            if config["PE_MODERN"]:
                target_list.append(f"fastq_inputs/PE_mod/{config['sample_prefix']}_R1_adRm.fastq.gz")
            elif config["PE_ANCIENT"]:
                target_list.append(f"fastq_inputs/PE_anc/{config['sample_prefix']}_adRm.fastq.gz")
            elif config["SE"]:
                target_list.append(f"fastq_inputs/SE/{config['sample_prefix']}_adRm.fastq.gz")

        config_sample = os.path.join(args.sample_output_dir, "database_fetch_config.yaml")

        print(config)
        with open(config_sample, "w") as fout:
            yaml.safe_dump(config, fout, default_flow_style=False)

        config["workflow_dir"] = os.path.join(BASE_DIR, "workflow")  # TODO tidy up

        target_list = [os.path.join(args.sample_output_dir, target) for target in target_list]
        snakefile = os.path.join(BASE_DIR, "workflow/sample.smk")

        return self._run_snakemake(snakefile, args, config, target_list)




    def analyse(self):
        parser = argparse.ArgumentParser(description="Analyse a sample")

        parser.add_argument(
            "-m",
            "--mode",
            choices=["filter", "align", "likelihoods", "probabilities", "abundances", "reads", "mapdamage",],
            help="Analysis mode for the selected sample",
            metavar="",
        )
        parser.add_argument(
            "-D", "--database", help="Path to the database output directory. MANDATORY", metavar="",
        )
        parser.add_argument(
            "-S", "--sample", help="Path to the sample output directory. MANDATORY", metavar="",
        )
        parser.add_argument(
            "-g",
            "--genera",
            nargs="+",
            help="List containing the names of specific genera "
            "the abundances should be calculated "
            "on, separated by a space character <genus1 genus2 genus3 ...>",
            metavar="",
            default=[],
        )
        parser.add_argument(
            "-o", "--output", help="Path to results directory.", metavar="", dest="analysis_output_dir",
        )
        parser.add_argument(
            "-T",
            "--read-probability-threshold",
            help="Posterior probability threshold for a read to belong to a certain species. "
            "Chose from 0.5, 0.75 and 0.95 (default:0.75).",
            choices=[0.5, 0.75, 0.95],
            default=float(0.75),
            type=float,
            metavar="",
        )
        parser.add_argument(
            "-c", "--cores", help="Number of cores for HAYSTACK to use", metavar="", type=int, default=MAX_CPU,
        )
        parser.add_argument(
            "-M",
            "--mem",
            help="Max memory resources allowed to be used ofr indexing the input for "
            "the filtering alignment "
            "(default: max available memory {})".format(MAX_MEM_MB),
            type=float,
            default=MAX_MEM_MB,
            metavar="",
        )
        parser.add_argument(
            "-u",
            "--unlock",
            action="store_true",
            help="Unlock the working directory after smk is " "abruptly killed  <bool> (default: False)",
        )
        parser.add_argument(
            "-d", "--debug", action="store_true", help="Debug the HAYSTACK workflow <bool> (default: False)",
        )
        parser.add_argument(
            "-smk", "--snakemake", help="Snakemake flags (default: '')", metavar="",
        )
        parser.add_argument("--dry-run", action="store_true")

        argcomplete.autocomplete(parser)
        args = parser.parse_args(sys.argv[2:])

        if len(sys.argv) == 2:
            parser.print_help()
            parser.exit()

        print("The selected mode for sample analysis is {}".format(args.mode))

        snakefile = os.path.join(BASE_DIR, "workflow", "analyse.smk")
        if not os.path.exists(snakefile):
            sys.stderr.write("Error: cannot find Snakefile at {}\n".format(snakefile))
            sys.exit(-1)

        repo_config_file = os.path.join(BASE_DIR, "config", "config.yaml")
        user_config_file = os.path.join(str(Path.home()), ".haystack", "config.yaml")

        if not os.path.exists(user_config_file):
            raise ValidationError(
                "Please run haystack config first in order to set up your "
                "email address and desired path for storing the downloaded genomes."
            )

        with open(repo_config_file) as fin:
            self.config_default = yaml.safe_load(fin)

        with open(user_config_file) as fin:
            self.config_user = yaml.safe_load(fin)

        self.config_default.update((k, v) for k, v in self.config_user.items())

        if not os.path.exists(args.database):
            raise ValidationError(
                "The path you provided for the database output directory is not valid. " "Please provide a valid path."
            )

        if not os.path.exists(args.sample):
            raise ValidationError(
                "The path you provided for the sample related output is not valid. " "Please provide a valid path."
            )

        db_fetch_yaml = os.path.join(args.database, "database_fetch_config.yaml")
        if os.path.exists(db_fetch_yaml):
            with open(db_fetch_yaml) as fin:
                database_config = yaml.safe_load(fin)

        db_build_yaml = os.path.join(args.database, "database_build_config.yaml")
        if os.path.exists(db_build_yaml):
            with open(db_build_yaml) as fin:
                database_config = yaml.safe_load(fin)

        if os.path.exists(db_build_yaml) and os.path.exists(db_fetch_yaml):
            raise ValidationError("The database has not been build correctly. Please re build the database.")

        if os.path.exists(db_build_yaml) is False and os.path.exists(db_fetch_yaml) is False:
            raise ValidationError("The database has not been build correctly or at all. Please re build the database.")

        sample_yaml = os.path.join(args.sample, "sample_config.yaml")
        if os.path.exists(sample_yaml):
            with open(sample_yaml) as fin:
                sample_config = yaml.safe_load(fin)
        else:
            raise ValidationError(
                "The sample yaml file does not exist in the path you provided. "
                "Please make sure you have provided the right sample output path, or "
                "make sure that you have run haystack sample first. "
            )

        analysis_args = vars(args)

        analysis_config = {k: v for k, v in self.config_default.items()}
        analysis_config.update((k, v) for k, v in database_config.items())
        analysis_config.update((k, v) for k, v in sample_config.items())
        analysis_config.update((k, v) for k, v in analysis_args.items())

        analysis_config["analysis_output_dir"] = os.path.abspath(analysis_config["analysis_output_dir"])

        if analysis_config["analysis_output_dir"]:
            if os.path.exists(analysis_config["analysis_output_dir"]):
                if not os.access(analysis_config["analysis_output_dir"], os.W_OK):
                    raise ValidationError(
                        "This directory path you have provided is not writable. "
                        "Please chose another path for your sample output directory."
                    )
            else:
                if not os.access(os.path.dirname(analysis_config["analysis_output_dir"]), os.W_OK):
                    raise ValidationError(
                        "This directory path you have provided is not writable. "
                        "Please chose another path for your sample output directory."
                    )

        if analysis_config["analysis_output_dir"] is None:
            raise ValidationError(
                "Please provide a valid directory path for the species identification related outputs. "
                "If the directory does not exist, do not worry the method will create it."
            )
        # print(analysis_config)

        target_list = []

        if args.mode == "filter":
            bowtie = ""
            if analysis_config["PE_MODERN"]:
                bowtie = analysis_config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq_pair.readlen".format(
                    sample=analysis_config["sample_prefix"]
                )
            elif analysis_config["PE_ANCIENT"] or analysis_config["SE"]:
                bowtie = analysis_config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.readlen".format(
                    sample=analysis_config["sample_prefix"]
                )
            target_list.append(bowtie)

        if args.mode == "align":
            target_list.append(
                analysis_config["analysis_output_dir"]
                + "/sigma/{sample}_alignments.done".format(sample=analysis_config["sample_prefix"])
            )

        if args.mode == "likelihoods":
            target_list.append(
                analysis_config["analysis_output_dir"]
                + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv".format(
                    sample=analysis_config["sample_prefix"]
                )
            )

        if args.mode == "probabilities":
            target_list.append(
                analysis_config["analysis_output_dir"]
                + "/probabilities/{sample}/{sample}_posterior_probabilities.tsv".format(
                    sample=analysis_config["sample_prefix"]
                )
            )

        if args.mode == "abundances":
            target_list.append(
                analysis_config["analysis_output_dir"]
                + "/probabilities/{sample}/{sample}_posterior_abundance.tsv".format(
                    sample=analysis_config["sample_prefix"]
                )
            )

        if args.mode == "reads":
            target_list.append(
                analysis_config["analysis_output_dir"]
                + "/dirichlet_reads/{sample}_dirichlet_reads.done".format(sample=analysis_config["sample_prefix"])
            )

        if args.mode == "mapdamage":
            target_list.append(
                analysis_config["analysis_output_dir"]
                + "/mapdamage/{sample}_mapdamage.done".format(sample=analysis_config["sample_prefix"])
            )

        analysis_yaml = os.path.join(
            str(Path.home()), analysis_config["analysis_output_dir"], analysis_config["sample_prefix"] + "_config.yaml",
        )

        if not os.path.exists(analysis_yaml):
            os.makedirs(
                os.path.join(str(Path.home()), analysis_config["analysis_output_dir"]), exist_ok=True,
            )
            analysis_options = {k: v for k, v in analysis_config.items() if (k, v) not in self.config_default.items()}
            with open(analysis_yaml, "w") as outfile:
                yaml.safe_dump(analysis_options, outfile, default_flow_style=False)

        analysis_config["workflow_dir"] = os.path.join(BASE_DIR, "workflow")

        user_options = {k: v for k, v in analysis_args.items() if (k, v) not in self.config_default.items()}

        print("--------")
        print("RUN DETAILS")
        print("\n\tanalyse.smk: {}".format(snakefile))
        print("\n\tConfig Parameters:\n")
        if args.debug:
            for (key, value,) in analysis_config.items():
                print(f"{key:35}{value}")
        else:
            for (key, value,) in user_options.items():
                print(f"{key:35}{value}")

        print("\n\tTarget Output Files:\n")
        for target in target_list:
            print(target)
        print("--------")

        if args.debug:
            printshellcmds = True
            keepgoing = False
            restart_times = 0
        else:
            printshellcmds = False
            keepgoing = True
            restart_times = 3

        status = snakemake.snakemake(
            snakefile,
            config=analysis_config,
            targets=target_list,
            printshellcmds=printshellcmds,
            dryrun=args.dry_run,
            cores=int(args.cores),
            keepgoing=keepgoing,
            restart_times=restart_times,
            unlock=args.unlock,
            show_failed_logs=args.debug,
            resources={"entrez_api": MAX_ENTREZ_REQUESTS},
            use_conda=analysis_config["use_conda"],
        )

        # translate "success" into shell exit code of 0
        return 0 if status else 1

    @staticmethod
    def _run_snakemake(snakefile, args, config, target_list):
        """
        Helper function for running the snakemake workflow
        """
        print("HAYSTACK\n")
        print(f"Date: {datetime.datetime.now()}\n")

        print("Config parameters:\n")
        params = config if args.debug else vars(args)

        for key, value in params.items():
            if value or args.debug:
                print(f" {key}: {value}")
        print("\n")

        if args.debug:
            print("Target files:\n")
            for target in target_list:
                print(" " + target)
            print("\n")

        # get any extra snakemake params
        smk_params = config.pop("snakemake") or {}

        status = snakemake.snakemake(
            snakefile,
            config=config,
            targets=target_list,
            printshellcmds=args.debug,
            cores=int(args.cores),
            keepgoing=(not args.debug),
            restart_times=0 if args.debug else RESTART_TIMES,
            unlock=args.unlock,
            show_failed_logs=args.debug,
            resources={"entrez_api": MAX_ENTREZ_REQUESTS},
            use_conda=config["use_conda"],
            **smk_params,
        )

        # translate "success" into shell exit code of 0
        return 0 if status else 1


if __name__ == "__main__":
    Haystack()
