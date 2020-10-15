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
    JsonType
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


def interactive_config_input():
    entrez_email = ""
    count = 0
    while count < 3:
        entrez_email = input(
            "Please enter a valid email address. "
            "It is required, in order to access NCBI's Entrez API."
            "The address is stored locally on your computer only: "
        )
        if not re.match(r"[^@]+@[^@]+\.[^@]+", entrez_email):
            print("The email address you provided is not valid. Please try again")
            entrez_email = input(
                "Please enter a valid email address. "
                "It is required, in order to access NCBI's Entrez API. "
                "The address is stored locally on your computer only: "
            )
            count += 1
        else:
            break
    if count == 3:
        raise ValidationError("Please try haystack config again to input a valid email address.")

    cache = input(
        "Enter your preferred path for the genome cache. " "Press enter if you'd like to use the default location: "
    ) or os.path.join(
        str(Path.home()), "haystack/cache/"
    )  # TOOD use the fucking default!!

    batchsize = input(
        "Enter your preferred batchsize for fetching accession data from the NCBI. "
        "Press enter if you'd like to use the default value: "
    ) or int(5)

    mismatch_probability = input(
        "Enter your preferred mismatch probability. " "Press enter if you'd like to use the default value: "
    ) or float(0.05)

    bowtie2_threads = input(
        "Enter your preferred number of threads that bowtie2 can use. "
        "Press enter if you'd like to use the default value: "
    ) or int(5)

    bowtie2_scaling = input(
        "Enter your preferred scaling factor for the size of the bowtie2 index chunks. "
        "Press enter if you'd like to use the default value: "
    ) or float(15)

    use_conda = (
        input(
            "Enter your preference about using conda as a package manager. "
            "Press enter if you'd like to use the default value: "
        )
        or True
    )

    cache = cache.rstrip("/")

    if os.path.exists(cache):
        if not os.access(cache, os.W_OK):
            raise ValidationError(
                "This directory path you have provided is not writable. "
                "Please chose another path for your genomes directory."
            )
    else:
        if not os.access(os.path.dirname(cache), os.W_OK):
            raise ValidationError(
                "This directory path you have provided is not writable. "
                "Please chose another path for your genomes directory."
            )

    user_data = {
        "cache": os.path.abspath(cache),
        "email": entrez_email,
        "batchsize": int(batchsize),
        "mismatch_probability": float(mismatch_probability),
        "bowtie2_threads": int(bowtie2_threads),
        "bowtie2_scaling": float(bowtie2_scaling),
        "use_conda": use_conda,
    }

    check_config_arguments(user_data)

    return user_data


def check_config_arguments(args):
    """Function to check config arguments and raise errors if they are not suitable"""

    if args["email"]:
        if not re.match(r"[^@]+@[^@]+\.[^@]+", args["email"]):
            print("The email address you provided is not valid. Please try again")
            args["email"] = input(
                "Please enter a valid email address. "
                "It is required, in order to access NCBI's Entrez API. "
                "The address is stored locally on your computer only: "
            )

    if args["cache"]:
        if os.path.exists(args["cache"]):
            if not os.access(args["cache"], os.W_OK):
                raise ValidationError(
                    "This directory path you have provided is not writable. "
                    "Please chose another path for your genomes directory."
                )
        else:
            if not os.access(os.path.dirname(args["cache"]), os.W_OK):
                raise ValidationError(
                    "This directory path you have provided is not writable. "
                    "Please chose another path for your genomes directory."
                )

    if args["batchsize"]:
        if not isinstance(args["batchsize"], int):
            raise ValidationError("Please provide a positive integer for batchsize.")
        if not args["batchsize"] > 0:
            raise ValidationError("Please provide a positive integer for batchsize.")

    if args["mismatch_probability"]:
        if not (args["mismatch_probability"], float):
            raise ValidationError("Please provide a positive float for mismatch probability.")
        if not args["mismatch_probability"] > 0:
            raise ValidationError("Please provide a positive float for mismatch probability.")

    if args["bowtie2_threads"]:
        if not isinstance(args["bowtie2_threads"], int):
            raise ValidationError("Please provide a positive integer for the bowtie2 threads.")
        if not args["bowtie2_threads"] > 0:
            raise ValidationError("Please provide a positive integer for the bowtie2 threads.")

    if args["bowtie2_scaling"]:
        if not (args["bowtie2_scaling"], float):
            raise ValidationError("Please provide a positive float for the bowtie2 scaling factor.")
        if not args["bowtie2_scaling"] > 0:
            raise ValidationError("Please provide a positive float for the bowtie2 scaling factor.")

    if args["use_conda"]:
        if args["use_conda"] not in ["True", "False"]:
            raise ValidationError("Please either specify True or False for using conda.")


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("True", "true"):
        return True
    elif v.lower() in ("False", "false"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------


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
        )

        parser.add_argument(
            "--batchsize",
            help=f"Batch size for fetching records from NCBI (default: {self.config_default['batchsize']})",
            type=PositiveIntType(),
            metavar="<int>",
        )

        parser.add_argument(
            "--mismatch-probability",
            help=f"Base mismatch probability (default: {self.config_default['mismatch_probability']})",
            type=FloatRangeType(0.01, 0.10),
            metavar="<float>",
        )

        parser.add_argument(
            "--bowtie2-threads",
            help=f"Number of threads to use for each bowtie2 alignment "
            f"(default: {self.config_default['bowtie2_threads']})",
            type=IntRangeType(1, MAX_CPU),
            metavar="<int>",
        )

        parser.add_argument(
            "--bowtie2-scaling",
            help=f"Rescaling factor to keep the bowtie2 mutlifasta index below the maximum memory limit "
            f"(default: {self.config_default['bowtie2_scaling']})",
            type=FloatRangeType(0, 100),
            metavar="<float>",
        )

        parser.add_argument(
            "--use-conda",
            help=f"Use conda as a package manger (default: {self.config_default['use_conda']})",
            type=BoolType(),
            metavar="<bool>",
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
        Add the common arguments shared by all commands
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
            help="Unlock the working directory following a crash or hard restart (default: False).",
            action="store_true",
        )

        parser.add_argument(
            "--debug", help="Enable debugging mode (default: False)", action="store_true",
        )

        parser.add_argument(
            "--snakemake",
            help="Pass additional flags to the `snakemake` scheduler.",
            metavar="'<json>'",
            type=JsonType()
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
        if (
            args.refseq_rep is False
            and args.query is None
            and args.query_file is None
            and args.accessions is None
            and args.sequences is None
        ):
            raise ValidationError(
                "Please specify at least one of --refseq-rep, --query, --query-file, --accessions or --sequences"
            )

        if args.mtDNA and args.refseq_rep:
            raise ValidationError(
                "Please specify either `--mtDNA` or `--refseq-rep` but not both."
            )

        if args.query and args.query_file:
            raise ValidationError(
                "Please specify either `--query <query>` or `--query-file <path>` but not both."
            )

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

        config['workflow_dir'] = os.path.join(BASE_DIR, "workflow")  # TODO tidy up
        config['mtDNA'] = str(args.mtDNA).lower()
        config["sequences"] = config["sequences"] or ""

        target_list = [os.path.join(args.db_output, target) for target in target_list]
        snakefile = os.path.join(BASE_DIR, "workflow/database.smk")

        return self._run_snakemake(snakefile, args, config, target_list)

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
        smk_params = config.pop('snakemake', {})

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
            use_conda=config['use_conda'],
            **smk_params
        )

        # translate "success" into shell exit code of 0
        return 0 if status else 1

    def sample(self):
        parser = argparse.ArgumentParser(description="Prepare a sample for analysis")

        parser.add_argument(
            "-p",
            "--sample-prefix",
            help="Sample prefix for all the future analysis. Optional if SRA accession is provided instead" " <str>",
            metavar="",
            default="",
        )
        parser.add_argument(
            "-o",
            "--output",
            help="Path to the directory where all the sample related outputs are going to be stored <str>",
            metavar="",
            default="",
            dest="sample_output_dir",
        )

        parser.add_argument(
            "-f", "--fastq", help="Path to the fastq input file. Can be raw or with adapters removed", metavar="",
        )
        parser.add_argument(
            "-f1",
            "--fastq-r1",
            help="Path to the mate 1 fastq input file, if reads are PE. " "Can be raw or with adapters removed",
            metavar="",
        )
        parser.add_argument(
            "-f2",
            "--fastq-r2",
            help="Path to the mate 2 fastq input file, if reads are PE. " "Can be raw or with adapters removed",
            metavar="",
        )

        parser.add_argument(
            "-SA",
            "--sra",
            help="Fetch raw data files from the SRA using the provided accession code <str>",
            metavar="",
        )

        parser.add_argument(
            "-C",
            "--collapse",
            help="Collapse paired end reads <bool> (default: False)",
            default=False,
            action="store_true",
        )
        parser.add_argument(
            "-T",
            "--not-trim-adapters",
            help="Do not remove adapters from raw fastq files <bool> (default: False)",
            action="store_true",
        )
        parser.add_argument(
            "-TF",
            "--adaperremoval-flags",
            help="Additional flags to provide to Adapterremoval <str>",
            default="",
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

        snakefile = os.path.join(BASE_DIR, "workflow", "sample.smk")
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

        sample_args = vars(args)
        # print(sample_args)
        sample_config = {k: v for k, v in self.config_default.items()}
        sample_config.update((k, v) for k, v in sample_args.items())

        sample_config["sample_output_dir"] = os.path.abspath(sample_config["sample_output_dir"])

        if sample_config["sample_output_dir"]:
            if os.path.exists(sample_config["sample_output_dir"]):
                if not os.access(sample_config["sample_output_dir"], os.W_OK):
                    raise ValidationError(
                        "This directory path you have provided is not writable. "
                        "Please chose another path for your sample output directory."
                    )
            else:
                if not os.access(os.path.dirname(sample_config["sample_output_dir"]), os.W_OK) and not os.access(
                    os.path.dirname(os.path.dirname(sample_config["sample_output_dir"])), os.W_OK,
                ):
                    raise ValidationError(
                        "This directory path you have provided is not writable. "
                        "Please chose another path for your sample output directory."
                    )

        if sample_config["sample_output_dir"] is None:
            raise ValidationError(
                "Please provide a valid directory path for the sample related outputs. "
                "If the directory does not exist, do not worry the method will create it."
            )

        sample_config["PE_ANCIENT"] = False
        sample_config["PE_MODERN"] = False
        sample_config["SE"] = False

        if sample_config["fastq_r1"] is not None and sample_config["fastq_r2"] is not None:
            if sample_config["collapse"]:
                sample_config["PE_ANCIENT"] = True
            else:
                sample_config["PE_MODERN"] = True

        if sample_config["fastq"] is not None:
            sample_config["SE"] = True

        if (sample_config["fastq_r1"] or sample_config["fastq_r2"] is not None) and (
            sample_config["fastq"] is not None
        ):
            raise ValidationError("Please use a correct combination of PE or SE reads.")

        if sample_config["fastq"]:
            if not os.path.exists(sample_config["fastq"]):
                raise ValidationError(
                    "The file path you provided to --fastq does not exist. " "Please provide a valid path"
                )

        if sample_config["fastq_r1"]:
            if not os.path.exists(sample_config["fastq_r1"]):
                raise ValidationError(
                    "The file path you provided to --fastq-r1 does not exist. " "Please provide a valid path"
                )

        if sample_config["fastq_r2"]:
            if not os.path.exists(sample_config["fastq_r2"]):
                raise ValidationError(
                    "The file path you provided to --fastq-r2 does not exist. " "Please provide a valid path"
                )

        if sample_config["sra"] is None and sample_config["sample_prefix"] is None:
            raise ValidationError("Please provide a prefix name for the sample you want to analyse.")

        if sample_config["sra"] is not None:
            sample_config["sample_prefix"] = sample_config["sra"]

            Entrez.email = sample_config["email"]
            sra_id = Entrez.read(Entrez.esearch(db="sra", term=sample_config["sra"]))["IdList"]
            if "paired" in str(Entrez.read(Entrez.esummary(db="sra", id=sra_id))).lower():
                if sample_config["collapse"]:
                    sample_config["PE_ANCIENT"] = True
                else:
                    sample_config["PE_MODERN"] = True
            else:
                sample_config["SE"] = True

        if sample_config["sra"] is not None:
            if sample_config["PE_MODERN"] or sample_config["PE_ANCIENT"]:
                sample_config["fastq_r1"] = "{prefix}/sra_data/PE/{accession}_R1.fastq.gz".format(
                    prefix=sample_config["sample_output_dir"], accession=sample_config["sra"],
                )
                sample_config["fastq_r2"] = "{prefix}/sra_data/PE/{accession}_R2.fastq.gz".format(
                    prefix=sample_config["sample_output_dir"], accession=sample_config["sra"],
                )
            elif sample_config["SE"]:
                sample_config["fastq"] = "{prefix}/sra_data/SE/{accession}.fastq.gz".format(
                    prefix=sample_config["sample_output_dir"], accession=sample_config["sra"],
                )

        if sample_config["fastq"] and sample_config["collapse"]:
            raise (
                ValidationError(
                    "You cannot collapse SE reads. Please delete the --collapse flag from your command, "
                    "or provide a different set of input files."
                )
            )

        target_list = [
            sample_config["sample_output_dir"]
            + "/fastq_inputs/meta/{sample}.size".format(sample=sample_config["sample_prefix"])
        ]

        if sample_config["not_trim_adapters"]:
            sample_config["trim_adapters"] = False
        else:
            sample_config["trim_adapters"] = True

        data_preprocessing = ""
        if sample_config["trim_adapters"]:
            if sample_config["PE_MODERN"]:
                data_preprocessing = sample_config[
                    "sample_output_dir"
                ] + "/fastq_inputs/PE_mod/{sample}_R1_adRm.fastq.gz".format(sample=sample_config["sample_prefix"])
            elif sample_config["PE_ANCIENT"]:
                data_preprocessing = sample_config[
                    "sample_output_dir"
                ] + "/fastq_inputs/PE_anc/{sample}_adRm.fastq.gz".format(sample=sample_config["sample_prefix"])
            elif sample_config["SE"]:
                data_preprocessing = sample_config[
                    "sample_output_dir"
                ] + "/fastq_inputs/SE/{sample}_adRm.fastq.gz".format(sample=sample_config["sample_prefix"])
            target_list.append(data_preprocessing)

        sample_yaml = os.path.join(str(Path.home()), sample_config["sample_output_dir"], "sample_config.yaml",)

        if not os.path.exists(sample_yaml):
            os.makedirs(
                os.path.join(str(Path.home()), sample_config["sample_output_dir"]), exist_ok=True,
            )
            sample_options = {k: v for k, v in sample_config.items() if (k, v) not in self.config_default.items()}
            with open(sample_yaml, "w") as outfile:
                yaml.safe_dump(sample_options, outfile, default_flow_style=False)

        sample_config["workflow_dir"] = os.path.join(BASE_DIR, "workflow")

        user_options = {k: v for k, v in sample_args.items() if (k, v) not in self.config_default.items()}
        # print(database_config)
        # print(sample_config)

        print("--------")
        print("RUN DETAILS")
        print("\n\tSnakefile: {}".format(snakefile))
        print("\n\tConfig Parameters:\n")
        if args.debug:
            for (key, value,) in sample_config.items():
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
            config=sample_config,
            targets=target_list,
            printshellcmds=printshellcmds,
            dryrun=args.dry_run,
            cores=int(args.cores),
            keepgoing=keepgoing,
            restart_times=restart_times,  # TODO find a better solution to this... 15 is way too many!
            unlock=args.unlock,
            show_failed_logs=args.debug,
            resources={"entrez_api": MAX_ENTREZ_REQUESTS},
            use_conda=sample_config["use_conda"],
        )

        # translate "success" into shell exit code of 0
        return 0 if status else 1

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


if __name__ == "__main__":
    Haystack()
