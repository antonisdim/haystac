#! /usr/bin/env python
"""
Execution script for snakemake workflows.
"""
__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

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
    EmailType,
    WritablePathType,
    PositiveIntType,
    FloatRangeType,
    IntRangeType,
    BoolType,
)

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

MEGABYTE = float(1024 ** 2)
MAX_MEM_MB = virtual_memory().total / MEGABYTE
MAX_CPU = cpu_count()

CONFIG_DEFAULT = "./config/config.yaml"
CONFIG_USER = pathlib.Path("~/.haystack/config.yaml").expanduser()

thisdir = os.path.abspath(os.path.dirname(__file__))

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
        raise RuntimeError("Please try haystack config again to input a valid email address.")

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
            raise RuntimeError(
                "This directory path you have provided is not writable. "
                "Please chose another path for your genomes directory."
            )
    else:
        if not os.access(os.path.dirname(cache), os.W_OK):
            raise RuntimeError(
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
                raise RuntimeError(
                    "This directory path you have provided is not writable. "
                    "Please chose another path for your genomes directory."
                )
        else:
            if not os.access(os.path.dirname(args["cache"]), os.W_OK):
                raise RuntimeError(
                    "This directory path you have provided is not writable. "
                    "Please chose another path for your genomes directory."
                )

    if args["batchsize"]:
        if not isinstance(args["batchsize"], int):
            raise RuntimeError("Please provide a positive integer for batchsize.")
        if not args["batchsize"] > 0:
            raise RuntimeError("Please provide a positive integer for batchsize.")

    if args["mismatch_probability"]:
        if not (args["mismatch_probability"], float):
            raise RuntimeError("Please provide a positive float for mismatch probability.")
        if not args["mismatch_probability"] > 0:
            raise RuntimeError("Please provide a positive float for mismatch probability.")

    if args["bowtie2_threads"]:
        if not isinstance(args["bowtie2_threads"], int):
            raise RuntimeError("Please provide a positive integer for the bowtie2 threads.")
        if not args["bowtie2_threads"] > 0:
            raise RuntimeError("Please provide a positive integer for the bowtie2 threads.")

    if args["bowtie2_scaling"]:
        if not (args["bowtie2_scaling"], float):
            raise RuntimeError("Please provide a positive float for the bowtie2 scaling factor.")
        if not args["bowtie2_scaling"] > 0:
            raise RuntimeError("Please provide a positive float for the bowtie2 scaling factor.")

    if args["use_conda"]:
        if args["use_conda"] not in ["True", "False"]:
            raise RuntimeError("Please either specify True or False for using conda.")


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
            self.config_default["cache"] = pathlib.Path(self.config_default["cache"]).expanduser()

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

        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def config(self):
        """
        Advanced configuration options for haystack
        """
        parser = argparse.ArgumentParser(description="Advanced configuration options for haystack")

        parser.add_argument(
            # "-e",
            "--email",
            help="Email address for NCBI identification (mandatory).",
            metavar="<address>",
            type=EmailType("RFC5322"),
        )

        parser.add_argument(
            # "-c",
            "--cache",
            help=f"Cache folder for all the genomes downloaded from NCBI (default: {self.config_default['cache']}).",
            metavar="<path>",
            type=WritablePathType(),
        )

        parser.add_argument(
            # "-b",
            "--batchsize",
            help=f"Batch size for fetching records from NCBI (default: {self.config_default['batchsize']})",
            type=PositiveIntType(),
            metavar="<int>",
        )

        parser.add_argument(
            # "-mp",
            "--mismatch-probability",
            help=f"Base mismatch probability (default: {self.config_default['mismatch_probability']})",
            type=FloatRangeType(0.01, 0.10),
            metavar="<float>",
        )

        parser.add_argument(
            # "-t",
            "--bowtie2-threads",
            help=f"Number of threads to use for each bowtie2 alignment "
            f"(default: {self.config_default['bowtie2_threads']})",
            type=IntRangeType(1, MAX_CPU),
            metavar="<int>",
        )

        parser.add_argument(
            # "-s",
            "--bowtie2-scaling",
            help=f"Rescaling factor to keep the bowtie2 mutlifasta index below the maximum memory limit "
            f"(default: {self.config_default['bowtie2_scaling']})",
            type=FloatRangeType(0, 100),
            metavar="<float>",
        )

        parser.add_argument(
            # "-cn",
            "--use-conda",
            help=f"Use conda as a package manger (default: {self.config_default['use_conda']})",
            type=BoolType(),
            metavar="<bool>",
        )

        # now that we're inside a subcommand, ignore the first two argvs
        argcomplete.autocomplete(parser)
        args = parser.parse_args(sys.argv[2:])

        config_user = dict()

        # get the user choices that differ from the defaults
        for key, value in vars(args).items():
            if value != self.config_default.get(key) and value is not None:
                config_user[key] = value

        # save the user config
        with open(CONFIG_USER, "w") as fout:
            yaml.safe_dump(config_user, fout)

    def database(self):
        parser = argparse.ArgumentParser(description="Build the database for haystack to use")

        parser.add_argument("--dry-run", action="store_true")

        parser.add_argument(
            "-m",
            "--mode",
            choices=["fetch", "index", "build"],
            help="Database creation mode for haystack",
            metavar="",
            default="build",
        )

        parser.add_argument(
            "-o", "--output", help="Path to the database output directory.", metavar="", dest="db_output",
        )
        parser.add_argument(
            "-R",
            "--refseq-rep",
            help="Use the prokaryotic representative species of the RefSeq DB "
            "for the species id pipeline. only species no strains. "
            "either or both of --refseq-rep and "
            "--query should be set (default: False)",
            type=bool,
            default=False,
            metavar="",
        )
        parser.add_argument(
            "-MT",
            "--mtDNA",
            help="Download mitochondrial genomes for eukaryotes only. "
            "Do not use with --refseq-rep or any queries for prokaryotes (default: False)",
            action="store_true",
        )
        parser.add_argument(
            "-q",
            "--query",
            help="Actual NCBI query in the NCBI query language. "
            "Please refer to the documentation on how to construct one correctly.",
            metavar="",
        )
        parser.add_argument(
            "-Q",
            "--query-file",  # todo check what Evan meant here
            help="Actual NCBI query in the NCBI query language, stored in a simple text file.",
            metavar="",
        )
        parser.add_argument(
            "-r",
            "--rank",
            help="Taxonomic rank to perform the identifications on (genus, species, subspecies, serotype) "
            "<str> (default: species)",
            choices=["genus", "species", "subspecies", "serotype"],
            default="species",
            metavar="",
        )
        parser.add_argument(
            "-s",
            "--sequences",
            help="TAB DELIMITED input file containing the the name of the taxon with no special characters, "
            "and an underscore '_' instead of spaces, a user defined accession code and the path of the fasta file. "
            "The fasta file that the path point to can be either uncompressed or compressed with gzip/bgzip",
            metavar="",
            default="",
        )

        parser.add_argument(
            "-a",
            "--accessions",
            help="TAB DELIMITED input file containing the the name of the taxon with no special characters, "
            "and an underscore '_' instead of spaces, a user defined valid NCBI nucleotide, assembly or WGS "
            "accession code. ",
            metavar="",
            default="",
        )

        parser.add_argument(
            "-S",
            "--seed",
            help="Seed for the randomization of the genomes that each index chunk will include <int> (default 1)",
            metavar="",
            type=int,
            default=int(1),
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
            "-c", "--cores", help="Number of cores for HAYSTACK to use", type=int, metavar="", default=MAX_CPU,
        )
        parser.add_argument(
            "-M",
            "--mem",
            help="Max memory resources allowed to be used for indexing the input for "
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
            "-smk", "--snakemake", help="Snakemake flags (default: '')", metavar="",  # todo don't know how to do that
        )

        # now that we're inside a subcommand, ignore the first
        # TWO argvs, ie the command (git) and the subcommand (commit)
        argcomplete.autocomplete(parser)
        args = parser.parse_args(sys.argv[2:])

        if len(sys.argv) == 2:
            parser.print_help()
            parser.exit()

        snakefile = os.path.join(thisdir, "workflow", "database.smk")
        if not os.path.exists(snakefile):
            sys.stderr.write("Error: cannot find Snakefile at {}\n".format(snakefile))
            sys.exit(-1)

        print("Running db creation, option={}".format(args.mode))

        # actual arg parsing

        repo_config_file = os.path.join(thisdir, "config", "config.yaml")
        user_config_file = os.path.join(str(Path.home()), ".haystack", "config.yaml")

        if not os.path.exists(user_config_file):
            raise RuntimeError(
                "Please run haystack config first in order to set up your "
                "email address and desired path for storing the downloaded genomes."
            )

        # rip_config = self.config
        with open(repo_config_file) as fin:
            repo_rip_config = yaml.safe_load(fin)

        with open(user_config_file) as fin:
            user_rip_config = yaml.safe_load(fin)

        repo_rip_config.update((k, v) for k, v in user_rip_config.items())

        database_args = vars(args)

        database_config = {k: v for k, v in repo_rip_config.items()}
        database_config.update((k, v) for k, v in database_args.items())

        if (
            database_config["refseq_rep"] is False
            and database_config["query"] is None
            and database_config["query_file"] is None
            and database_config["accessions"] is None
            and database_config["sequences"] is None
        ):
            raise RuntimeError(
                "Please specify where HAYSTACK should get the database sequences from "
                "(query, RefSeq Rep, custom accessions or custom seqeunces"
            )

        if database_config["mtDNA"] and database_config["refseq_rep"]:
            raise RuntimeError(
                "These flags are mutually exclusive. Either pick refseq-rep for prokaryotic "
                "related queries or pick mtDNA for eukaryotic related queries."
            )

        if database_config["query"] and database_config["query_file"]:
            raise RuntimeError(
                "You cannot provide both a query and a query file. "
                "Please chose only one option, and provide the only its respective flag."
            )

        if database_config["query_file"]:
            with open(database_config["query_file"], "r") as fin:
                database_config["query"] = fin.read().rstrip()

            if database_config["query"] == "" or database_config["query"] == " ":
                raise RuntimeError("The query file you provided was empty. Please provide a file with a valid query.")

        database_config["db_output"] = os.path.abspath(database_config["db_output"])

        if database_config["db_output"]:
            if os.path.exists(database_config["db_output"]):
                if not os.access(database_config["db_output"], os.W_OK):
                    raise RuntimeError(
                        "This directory path you have provided is not writable. "
                        "Please chose another path for your database output directory."
                    )
            else:
                if not os.access(os.path.dirname(database_config["db_output"]), os.W_OK):
                    raise RuntimeError(
                        "This directory path you have provided is not writable. "
                        "Please chose another path for your database output directory.."
                    )

        if database_config["db_output"] is None:
            raise RuntimeError(
                "Please provide a valid directory path for the database outputs. "
                "If the directory does not exist, do not worry the method will create it."
            )

        target_list = []

        if database_config["mode"] == "fetch":
            if database_config["query"] != "":
                target_list.append(database_config["db_output"] + "/bowtie/entrez_query.fasta.gz")
            if database_config["refseq_rep"]:
                target_list.append(database_config["db_output"] + "/bowtie/refseq_prok.fasta.gz")
            if database_config["sequences"] != "":
                target_list.append(database_config["db_output"] + "/bowtie/custom_seqs.fasta.gz")
            if database_config["accessions"] != "":
                target_list.append(database_config["db_output"] + "/bowtie/custom_acc.fasta.gz")

            database_fetch_yaml = os.path.join(
                str(Path.home()), database_config["db_output"], "database_fetch_config.yaml",
            )
            if not os.path.exists(database_fetch_yaml):
                os.makedirs(
                    os.path.join(str(Path.home()), database_config["db_output"]), exist_ok=True,
                )
                with open(database_fetch_yaml, "w") as outfile:
                    yaml.safe_dump(database_config, outfile, default_flow_style=False)

            print("Please run haystack database --mode index after this step.")

        if database_config["mode"] == "index":
            target_list.append(database_config["db_output"] + "/bowtie/bowtie_index.done")

            database_fetch_yaml = os.path.join(
                str(Path.home()), database_config["db_output"], "database_fetch_config.yaml",
            )
            if not os.path.exists(database_fetch_yaml):
                raise RuntimeError(
                    "Please run haystack database --mode fetch first, and then proceed indexing the database."
                )

            with open(database_fetch_yaml, "r") as fin:
                database_config = yaml.safe_load(fin)

        if database_config["mode"] == "build":
            target_list.append(database_config["db_output"] + "/idx_database.done")
            target_list.append(database_config["db_output"] + "/bowtie/bowtie_index.done")

            database_fetch_yaml = os.path.join(
                str(Path.home()), database_config["db_output"], "database_fetch_config.yaml",
            )
            if os.path.exists(database_fetch_yaml):
                raise RuntimeError(
                    "You can not run haystack database --mode build after running --mode fetch. "
                    "You need to run haystack database --mode index instead."
                )

            database_build_yaml = os.path.join(
                str(Path.home()), database_config["db_output"], "database_build_config.yaml",
            )

            if not os.path.exists(database_build_yaml):
                os.makedirs(
                    os.path.join(str(Path.home()), database_config["db_output"]), exist_ok=True,
                )
                with open(database_build_yaml, "w") as outfile:
                    yaml.safe_dump(database_config, outfile, default_flow_style=False)

        database_config["workflow_dir"] = os.path.join(thisdir, "workflow")
        database_config["mtDNA"] = str(database_config["mtDNA"]).lower()

        user_options = {k: v for k, v in database_args.items() if (k, v) not in repo_rip_config.items()}
        # print(database_config)

        print("--------")
        print("RUN DETAILS")
        print("\n\tSnakefile: {}".format(snakefile))
        print("\n\tConfig Parameters:\n")
        if database_config["debug"]:
            for (key, value,) in database_config.items():
                print(f"{key:35}{value}")
        else:
            for (key, value,) in user_options.items():
                print(f"{key:35}{value}")

        print("\n\tTarget Output Files:\n")
        for target in target_list:
            print(target)
        print("--------")

        if database_config["debug"]:
            printshellcmds = True
            keepgoing = False
            restart_times = 0
        else:
            printshellcmds = False
            keepgoing = True
            restart_times = 3

        status = snakemake.snakemake(
            snakefile,
            config=database_config,
            targets=target_list,
            printshellcmds=printshellcmds,
            dryrun=args.dry_run,
            cores=int(args.cores),
            keepgoing=keepgoing,
            restart_times=restart_times,  # TODO find a better solution to this... 15 is way too many!
            unlock=args.unlock,
            show_failed_logs=args.debug,
            resources={"entrez_api": MAX_ENTREZ_REQUESTS},
            use_conda=database_config["use_conda"],
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

        snakefile = os.path.join(thisdir, "workflow", "sample.smk")
        if not os.path.exists(snakefile):
            sys.stderr.write("Error: cannot find Snakefile at {}\n".format(snakefile))
            sys.exit(-1)

        repo_config_file = os.path.join(thisdir, "config", "config.yaml")
        user_config_file = os.path.join(str(Path.home()), ".haystack", "config.yaml")

        if not os.path.exists(user_config_file):
            raise RuntimeError(
                "Please run haystack config first in order to set up your "
                "email address and desired path for storing the downloaded genomes."
            )

        with open(repo_config_file) as fin:
            repo_rip_config = yaml.safe_load(fin)

        with open(user_config_file) as fin:
            user_rip_config = yaml.safe_load(fin)

        repo_rip_config.update((k, v) for k, v in user_rip_config.items())

        sample_args = vars(args)
        # print(sample_args)
        sample_config = {k: v for k, v in repo_rip_config.items()}
        sample_config.update((k, v) for k, v in sample_args.items())

        sample_config["sample_output_dir"] = os.path.abspath(sample_config["sample_output_dir"])

        if sample_config["sample_output_dir"]:
            if os.path.exists(sample_config["sample_output_dir"]):
                if not os.access(sample_config["sample_output_dir"], os.W_OK):
                    raise RuntimeError(
                        "This directory path you have provided is not writable. "
                        "Please chose another path for your sample output directory."
                    )
            else:
                if not os.access(os.path.dirname(sample_config["sample_output_dir"]), os.W_OK) and not os.access(
                    os.path.dirname(os.path.dirname(sample_config["sample_output_dir"])), os.W_OK,
                ):
                    raise RuntimeError(
                        "This directory path you have provided is not writable. "
                        "Please chose another path for your sample output directory."
                    )

        if sample_config["sample_output_dir"] is None:
            raise RuntimeError(
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
            raise RuntimeError("Please use a correct combination of PE or SE reads.")

        if sample_config["fastq"]:
            if not os.path.exists(sample_config["fastq"]):
                raise RuntimeError(
                    "The file path you provided to --fastq does not exist. " "Please provide a valid path"
                )

        if sample_config["fastq_r1"]:
            if not os.path.exists(sample_config["fastq_r1"]):
                raise RuntimeError(
                    "The file path you provided to --fastq-r1 does not exist. " "Please provide a valid path"
                )

        if sample_config["fastq_r2"]:
            if not os.path.exists(sample_config["fastq_r2"]):
                raise RuntimeError(
                    "The file path you provided to --fastq-r2 does not exist. " "Please provide a valid path"
                )

        if sample_config["sra"] is None and sample_config["sample_prefix"] is None:
            raise RuntimeError("Please provide a prefix name for the sample you want to analyse.")

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
                RuntimeError(
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
            sample_options = {k: v for k, v in sample_config.items() if (k, v) not in repo_rip_config.items()}
            with open(sample_yaml, "w") as outfile:
                yaml.safe_dump(sample_options, outfile, default_flow_style=False)

        sample_config["workflow_dir"] = os.path.join(thisdir, "workflow")

        user_options = {k: v for k, v in sample_args.items() if (k, v) not in repo_rip_config.items()}
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

        snakefile = os.path.join(thisdir, "workflow", "analyse.smk")
        if not os.path.exists(snakefile):
            sys.stderr.write("Error: cannot find Snakefile at {}\n".format(snakefile))
            sys.exit(-1)

        repo_config_file = os.path.join(thisdir, "config", "config.yaml")
        user_config_file = os.path.join(str(Path.home()), ".haystack", "config.yaml")

        if not os.path.exists(user_config_file):
            raise RuntimeError(
                "Please run haystack config first in order to set up your "
                "email address and desired path for storing the downloaded genomes."
            )

        with open(repo_config_file) as fin:
            repo_rip_config = yaml.safe_load(fin)

        with open(user_config_file) as fin:
            user_rip_config = yaml.safe_load(fin)

        repo_rip_config.update((k, v) for k, v in user_rip_config.items())

        if not os.path.exists(args.database):
            raise RuntimeError(
                "The path you provided for the database output directory is not valid. " "Please provide a valid path."
            )

        if not os.path.exists(args.sample):
            raise RuntimeError(
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
            raise RuntimeError("The database has not been build correctly. Please re build the database.")

        if os.path.exists(db_build_yaml) is False and os.path.exists(db_fetch_yaml) is False:
            raise RuntimeError("The database has not been build correctly or at all. Please re build the database.")

        sample_yaml = os.path.join(args.sample, "sample_config.yaml")
        if os.path.exists(sample_yaml):
            with open(sample_yaml) as fin:
                sample_config = yaml.safe_load(fin)
        else:
            raise RuntimeError(
                "The sample yaml file does not exist in the path you provided. "
                "Please make sure you have provided the right sample output path, or "
                "make sure that you have run haystack sample first. "
            )

        analysis_args = vars(args)

        analysis_config = {k: v for k, v in repo_rip_config.items()}
        analysis_config.update((k, v) for k, v in database_config.items())
        analysis_config.update((k, v) for k, v in sample_config.items())
        analysis_config.update((k, v) for k, v in analysis_args.items())

        analysis_config["analysis_output_dir"] = os.path.abspath(analysis_config["analysis_output_dir"])

        if analysis_config["analysis_output_dir"]:
            if os.path.exists(analysis_config["analysis_output_dir"]):
                if not os.access(analysis_config["analysis_output_dir"], os.W_OK):
                    raise RuntimeError(
                        "This directory path you have provided is not writable. "
                        "Please chose another path for your sample output directory."
                    )
            else:
                if not os.access(os.path.dirname(analysis_config["analysis_output_dir"]), os.W_OK):
                    raise RuntimeError(
                        "This directory path you have provided is not writable. "
                        "Please chose another path for your sample output directory."
                    )

        if analysis_config["analysis_output_dir"] is None:
            raise RuntimeError(
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
            analysis_options = {k: v for k, v in analysis_config.items() if (k, v) not in repo_rip_config.items()}
            with open(analysis_yaml, "w") as outfile:
                yaml.safe_dump(analysis_options, outfile, default_flow_style=False)

        analysis_config["workflow_dir"] = os.path.join(thisdir, "workflow")

        user_options = {k: v for k, v in analysis_args.items() if (k, v) not in repo_rip_config.items()}

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
