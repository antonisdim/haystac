#! /usr/bin/env python
"""
Execution script for snakemake workflows.
"""
__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import argparse
import datetime
import os
import shutil
import sys
from multiprocessing import cpu_count
from snakemake.dag import Batch

import snakemake
import yaml
from psutil import virtual_memory

from haystac import __version__
from haystac.workflow.scripts.entrez_utils import (
    ENTREZ_RATE_LOW,
    ENTREZ_RATE_HIGH,
)
from haystac.workflow.scripts.utilities import (
    ValidationError,
    ArgumentCustomFormatter,
    WritablePathType,
    FloatRangeType,
    IntRangeType,
    BoolType,
    JsonType,
    SequenceFileType,
    AccessionFileType,
    SraAccessionType,
    NuccoreQueryType,
    CheckExistingConfig,
    BatchType,
    FastqFile,
    PE,
    COLLAPSED,
    SE,
    md5,
    print_warning,
    print_error,
)

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

MAX_CPU = cpu_count()
MAX_MEM_MB = int(virtual_memory().total / 1024 ** 2)

CODE_DIR = os.path.abspath(os.path.dirname(__file__))
SNAKE_DIR = ".snakemake"

CONFIG_DEFAULT = f"{CODE_DIR}/config/config.yaml"
CONFIG_USER = os.path.abspath(os.path.expanduser("~/.haystac/config.yaml"))
CONFIG_RUNTIME = f"{SNAKE_DIR}/config.yaml"

COMMANDS = ["config", "database", "sample", "analyse"]

DATABASE_MODES = ["fetch", "index", "build"]
ANALYSIS_MODES = ["filter", "align", "likelihoods", "probabilities", "abundances", "reads"]
TAXONOMIC_RANKS = ["genus", "species", "subspecies", "serotype"]
REFSEQ_TABLES = ["prokaryote_rep", "viruses", "eukaryotes"]

# number to times to retry a rule that failed the first time
RESTART_TIMES = 2


class Haystac(object):
    """
    Command-line interface for running `haystac`
    """

    def __init__(self):
        parser = argparse.ArgumentParser(
            usage="""haystac <command> [<args>]

The haystac commands are:
   config         Configuration options
   database       Build a database of target species
   sample         Prepare a sample for analysis
   analyse        Analyse a sample against a database
   
""",
        )
        parser.add_argument("command", choices=COMMANDS, help="Command to run")

        # get the CLI arguments
        args = self._parse_args(parser, level=1)

        # load the config files
        self._load_config()

        try:
            # call the class method with the name given by `command`
            reval = getattr(self, args.command)()
            exit(reval)

        except ValidationError as error:
            print_error(f"{error}")

    @staticmethod
    def _parse_args(parser, level=1):
        """
        Handle the nested-subcommand CLI interface.
        """

        # if we don't have enough arguments, then print the help
        if len(sys.argv) <= level:
            parser.print_help()
            parser.exit()

        if {"-v", "--version"}.intersection(sys.argv):
            print(__version__)
            parser.exit()

        # slice the CLI arguments
        argv = sys.argv[1:2] if level == 1 else sys.argv[level:]

        return parser.parse_args(argv)

    def _load_config(self):
        """
        Load the default and user space config files.
        """

        # load the default config
        with open(CONFIG_DEFAULT) as fin:
            self.config_default = yaml.safe_load(fin)

            # resolve relative paths
            self.config_default["cache"] = os.path.abspath(os.path.expanduser(self.config_default["cache"]))

        try:
            # load the user config
            with open(CONFIG_USER) as fin:
                self.config_user = yaml.safe_load(fin)

        except FileNotFoundError:
            # make the config folder (if necessary)
            os.makedirs(os.path.dirname(CONFIG_USER), exist_ok=True)

            self.config_user = dict()

        # maximum concurrent Entrez requests
        self.max_entrez_requests = ENTREZ_RATE_HIGH if self.config_user.get("api_key") else ENTREZ_RATE_LOW

        # merge the config dictionaries (giving precedence to the config_user)
        self.config_merged = {**self.config_default, **self.config_user}

    def config(self):
        """
        Configuration options
        """
        # noinspection PyTypeChecker
        parser = argparse.ArgumentParser(
            prog="haystac config",
            description="Configuration options",
            formatter_class=ArgumentCustomFormatter,
            add_help=False,
        )

        optional = parser.add_argument_group("Optional arguments")

        # add the help option manually so we can control where it is shown in the menu
        optional.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        optional.add_argument("-v", "--version", action="help", help="Print the version number and exit")

        optional.add_argument(
            "--cache",
            help=f"Cache folder for storing genomes downloaded from NCBI and other shared data "
            f"(default: {self.config_default['cache']})",
            metavar="<path>",
            type=WritablePathType(),
            default=argparse.SUPPRESS,
        )

        optional.add_argument(
            "--clear-cache",
            help="Clear the contents of the cache folder, and delete the folder itself",
            action="store_true",
        )

        optional.add_argument(
            "--api-key",
            help=f"Personal NCBI API key (increases max concurrent requests from {ENTREZ_RATE_LOW} to "
            f"{ENTREZ_RATE_HIGH}, https://www.ncbi.nlm.nih.gov/account/register/)",
            metavar="<code>",
            default=argparse.SUPPRESS,
        )

        optional.add_argument(
            "--use-conda",
            help=f"Use conda as a package manger (default: {self.config_default['use_conda']})",
            type=BoolType(),
            metavar="<bool>",
            default=argparse.SUPPRESS,
        )

        # get the CLI arguments
        args = self._parse_args(parser, level=2)

        # get the user choices
        for key, value in vars(args).items():
            self.config_user[key] = value

        # resolve relative paths
        if self.config_user.get("cache"):
            self.config_user["cache"] = os.path.abspath(os.path.expanduser(self.config_user["cache"]))

        # clear the cache directory if specified
        if self.config_user.get("clear_cache"):
            shutil.rmtree(self.config_user["cache"])
            # delete the clear_cache key, so that the cache is not cleared every time unless the user specifies it
            del self.config_user["clear_cache"]

        # save the user config
        with open(CONFIG_USER, "w") as fout:
            yaml.safe_dump(self.config_user, fout, default_flow_style=False)

    def database(self):
        """
        Build a database of target species
        """
        # noinspection PyTypeChecker
        parser = argparse.ArgumentParser(
            prog="haystac database",
            description="Build a database of target species",
            formatter_class=ArgumentCustomFormatter,
            add_help=False,
        )

        required = parser.add_argument_group("Required arguments")

        required.add_argument(
            "--mode",
            help="Database creation mode for haystac [%(choices)s]",
            metavar="<mode>",
            choices=DATABASE_MODES,
            required=True,
        )

        required.add_argument(
            "--output",
            help="Path to the database output directory",
            metavar="<path>",
            dest="db_output",
            type=WritablePathType(),
            required=True,
        )

        choice = parser.add_argument_group("Required choice")

        choice.add_argument(
            "--query",
            help="Database query in the NCBI query language. "
            "Please refer to the documentation for assistance with constructing a valid query.",
            metavar="<query>",
            type=NuccoreQueryType(),
        )

        choice.add_argument(
            "--query-file",
            help="File containing a database query in the NCBI query language.",
            metavar="<path>",
            type=NuccoreQueryType(),
        )

        choice.add_argument(
            "--accessions-file",
            dest="accessions",
            help="Tab delimited file containing one record per row: the name of the taxon, "
            "and a valid NCBI accession code from the nucleotide, assembly or WGS databases.",
            metavar="<path>",
            type=AccessionFileType(),
        )

        choice.add_argument(
            "--sequences-file",
            dest="sequences",
            help="Tab delimited file containing one record per row: the name of the taxon, a user defined "
            "accession code, and the path to the fasta file (optionally compressed).",
            metavar="<path>",
            type=SequenceFileType(),
        )

        choice.add_argument(
            "--refseq-rep",
            help="Use one of the RefSeq curated tables to construct a DB. Includes all prokaryotic species "
            "(excluding strains) from the representative RefSeq DB, or all the species and strains from the "
            "viruses DB, or all the species and subspecies from the eukaryotes DB. "
            "If multiple accessions exist for a given species/strain, the first pair of species/accession is kept. "
            "Available RefSeq tables to use [%(choices)s].",
            choices=REFSEQ_TABLES,
            metavar="<table>",
            default=None,
        )

        optional = parser.add_argument_group("Optional arguments")

        optional.add_argument(
            "--force-accessions",
            help="Disable validation checks for 'anomalous' assembly flags in NCBI",
            action="store_true",
        )

        optional.add_argument(
            "--exclude-accessions",
            help="List of NCBI accessions to exclude.",
            metavar="<accession>",
            nargs="+",
            default=[],
        )

        optional.add_argument(
            "--resolve-accessions",
            help="Pick the first accession when two accessions for a taxon can be found in user provided input files",
            action="store_true",
        )

        optional.add_argument(
            "--bowtie2-scaling",
            help="Rescaling factor to keep the bowtie2 mutlifasta index below the maximum memory limit",
            type=FloatRangeType(0, 100),
            metavar="<float>",
            default=self.config_default["bowtie2_scaling"],
        )

        optional.add_argument(
            "--bowtie2-threads",
            help="Number of threads bowtie2 will use to index every individual genome in the database",
            type=IntRangeType(1, MAX_CPU),
            metavar="<int>",
            dest="bowtie2_threads_db",
            default=self.config_default["bowtie2_threads_db"],
        )

        optional.add_argument(
            "--rank",
            help="Taxonomic rank to perform the identifications on [%(choices)s]",
            metavar="<rank>",
            choices=TAXONOMIC_RANKS,
            default=self.config_default["rank"],
        )

        optional.add_argument(
            "--genera",
            help="List of genera to restrict the abundance calculations.",
            metavar="<genus>",
            nargs="+",
            default=None,
        )

        optional.add_argument(
            "--mtDNA",
            help="For eukaryotes, download mitochondrial genomes only. "
            "Not to be used with --refseq-rep or queries containing prokaryotes",
            action="store_true",
        )

        optional.add_argument(
            "--seed",
            help="Random seed for database indexing",
            metavar="<int>",
            type=int,
            default=self.config_default["seed"],
        )

        optional.add_argument(
            "--batch",
            help="Batch number for large`haystac database` workflows (e.g. --batch index_all_accessions=1/3). "
            "You will need to execute all batches before haystac is able to finish its workflow to the end.",
            metavar="<str>",
            type=BatchType(),
            default=None,
        )

        # add the common arguments
        self._common_arguments(parser)

        # get the CLI arguments
        args = self._parse_args(parser, level=2)

        # must specify at least one source for the database
        if args.mode != "index" and not (
            args.refseq_rep or args.query or args.query_file or args.accessions or args.sequences
        ):
            raise ValidationError(
                "Please specify at least one of --refseq-rep, --query, --query-file, --accessions or --sequences"
            )

        if args.mtDNA and args.refseq_rep:
            raise ValidationError("Please specify either `--mtDNA` or `--refseq-rep` but not both.")

        if args.query and args.query_file:
            raise ValidationError("Please specify either `--query <query>` or `--query-file <path>` but not both.")

        if args.query_file:
            # load the query file
            with open(args.query_file) as fin:
                args.query = fin.read().strip()

        # warn the user in case he is selecting only one thread for bt2 accession indexing
        if args.bowtie2_threads_db == 1:
            print_warning(
                "Please note that when bowtie2 uses only 1 thread for indexing its memory usage tends "
                "to be higher than when it is multithreaded. That means that HAYSTAC could end up using "
                "more memory than it is desired. Please be cautious when using only 1 thread "
                "for indexing individual genomes."
            )

        db_original = args.db_output

        # resolve relative paths
        args.db_output = os.path.abspath(os.path.expanduser(args.db_output))

        # cast all `None` paths as "" or else smk complains when parsing the rules
        args.query_file = args.query_file or ""
        args.accessions = args.accessions or ""
        args.sequences = args.sequences or ""

        # add all command line options to the merged config
        # noinspection PyDictCreation
        config = {**self.config_merged, **vars(args)}

        # store the md5 checksums for the database user input files
        config["query_file_md5"] = md5(args.query_file) if args.query_file else ""
        config["accessions_md5"] = md5(args.accessions) if args.query_file else ""
        config["sequences_md5"] = md5(args.sequences) if args.query_file else ""

        config_fetch = os.path.join(args.db_output, "database_fetch_config.yaml")
        config_build = os.path.join(args.db_output, "database_build_config.yaml")

        if args.refseq_rep:
            if not args.force_accessions:
                print_warning("Automatically setting `--force-accessions` in `--refseq-rep` mode")

            config["force_accessions"] = True

        target_list = list()

        if args.mode == "fetch":
            if args.query or args.accessions or args.sequences or args.refseq_rep or args.query_file:
                target_list.append("db_taxa_accessions.tsv")

            CheckExistingConfig(config_fetch, config)

            with open(config_fetch, "w") as fout:
                yaml.safe_dump(config, fout, default_flow_style=False)

        elif args.mode == "index":
            target_list.append("bowtie/bowtie_index.done")
            target_list.append("idx_database.done")

            try:
                with open(config_fetch, "r") as fin:
                    # config = {**yaml.safe_load(fin), **config}
                    fetch_args = {**yaml.safe_load(fin)}
                    fetch_args.update((k, v) for k, v in config.items() if v is not None)
                    config = fetch_args

                    CheckExistingConfig(config_fetch, config)
            except FileNotFoundError:
                raise ValidationError(
                    "Please run haystac `database --mode fetch` before attempting to index the database."
                )

        elif args.mode == "build":
            target_list.append("idx_database.done")
            target_list.append("bowtie/bowtie_index.done")
            target_list.append("db_taxa_accessions.tsv")

            if os.path.exists(config_fetch):
                raise ValidationError(
                    "Please run haystac `database --mode index` as the database has already been fetched."
                )

            if os.path.exists(config_build):
                CheckExistingConfig(config_build, config)

            with open(config_build, "w") as fout:
                yaml.safe_dump(config, fout, default_flow_style=False)

        target_list = [os.path.join(args.db_output, target) for target in target_list]

        # run the `haystac` workflow
        exit_code = self._run_snakemake("database", args, config, target_list)

        if args.mode == "fetch" and exit_code == 0 and not args.unlock:
            print(f"Please run `haystac database --mode index --output {db_original}` after this step.")

        return exit_code

    def sample(self):
        """
        Prepare a sample for analysis
        """
        # noinspection PyTypeChecker
        parser = argparse.ArgumentParser(
            prog="haystac sample",
            description="Prepare a sample for analysis",
            formatter_class=ArgumentCustomFormatter,
            add_help=False,
        )

        required = parser.add_argument_group("Required arguments")

        required.add_argument(
            "--output",
            help="Path to the sample output directory",
            metavar="<path>",
            dest="sample_output_dir",
            type=WritablePathType(),
            required=True,
        )

        choice = parser.add_argument_group("Required choice")

        choice.add_argument(
            "--fastq",
            help="Single-end fastq input file (optionally compressed).",
            metavar="<path>",
            type=FastqFile(),
        )

        choice.add_argument(
            "--fastq-r1",
            help="Paired-end forward strand (R1) fastq input file.",
            metavar="<path>",
            type=FastqFile(),
        )

        choice.add_argument(
            "--fastq-r2",
            help="Paired-end reverse strand (R2) fastq input file.",
            metavar="<path>",
            type=FastqFile(),
        )

        choice.add_argument(
            "--sra",
            help="Download fastq input from the SRA database",
            type=SraAccessionType(),
            metavar="<accession>",
        )

        optional = parser.add_argument_group("Optional arguments")

        optional.add_argument(
            "--collapse",
            help="Collapse overlapping paired-end reads, e.g. for aDNA",
            type=BoolType(),
            metavar="<bool>",
            default=self.config_default["collapse"],
        )

        optional.add_argument(
            "--trim-adapters",
            help="Automatically trim sequencing adapters from fastq input",
            type=BoolType(),
            metavar="<bool>",
            default=self.config_default["trim_adapters"],
        )

        # add the common arguments
        self._common_arguments(parser)

        # get the CLI arguments
        args = self._parse_args(parser, level=2)

        # must specify exactly one source for the sample
        if bool(args.fastq) + (bool(args.fastq_r1) and bool(args.fastq_r2)) + bool(args.sra) != 1:
            raise ValidationError("Please specify either --sra or --fastq, or both --fastq-r1 and --fastq-r2.")

        if (bool(args.fastq) and bool(args.fastq_r1)) or (bool(args.fastq) and bool(args.fastq_r2)):
            raise ValidationError(
                "Please specify either --fastq, or both --fastq-r1 and --fastq-r2. "
                "Any combination of these flags is not valid"
            )

        if bool(args.sra) and (bool(args.fastq_r1) or bool(args.fastq) or bool(args.fastq_r2)):
            raise ValidationError(
                "Please specify either --sra or --fastq, or both --fastq-r1 and --fastq-r2. "
                "Any combination of these flags is not valid"
            )

        if args.collapse and (not (args.fastq_r1 and args.fastq_r2) and (bool(args.sra) and args.sra[1] != "paired")):
            raise ValidationError("Collapse can only be used with --fastq-r1 and --fastq-r2 or sra PE data.")

        if args.collapse and not args.trim_adapters:
            raise ValidationError("Collapse can only be used with `--trim-adapters True`.")

        # resolve relative paths
        args.sample_output_dir = os.path.abspath(os.path.expanduser(args.sample_output_dir))

        # cast all `None` paths as "" or ele smk complains when parsing the rules
        args.fastq = args.fastq or ""
        args.fastq_r1 = args.fastq_r1 or ""
        args.fastq_r2 = args.fastq_r2 or ""

        # add all command line options to the merged config
        config = {**self.config_merged, **vars(args)}

        if not (args.fastq_r1 and args.fastq_r2):
            config["read_mode"] = SE
        if (args.fastq_r1 and args.fastq_r2) and args.collapse:
            config["read_mode"] = COLLAPSED
        if (args.fastq_r1 and args.fastq_r2) and not args.collapse:
            config["read_mode"] = PE

        # use the --output folder to name the sample
        config["sample_prefix"] = os.path.basename(config["sample_output_dir"].rstrip("/"))

        if args.sra:
            config["sra"], config["layout"] = args.sra

            # query the SRA to see if this is a paired-end library or not
            if config["layout"] == "paired":
                if config["collapse"]:
                    config["read_mode"] = COLLAPSED
                else:
                    config["read_mode"] = PE

                config["fastq_r1"] = config["sample_output_dir"] + f"/sra_data/PE/{config['sra']}_R1.fastq.gz"
                config["fastq_r2"] = config["sample_output_dir"] + f"/sra_data/PE/{config['sra']}_R1.fastq.gz"

            else:
                config["read_mode"] = SE
                config["fastq"] = config["sample_output_dir"] + f"/sra_data/SE/{config['sra']}.fastq.gz"

        target_list = list()
        target_list.append(f"fastq_inputs/meta/{config['sample_prefix']}.size")

        if config["trim_adapters"]:
            if config["read_mode"] == PE:
                target_list.append(f"fastq_inputs/{config['read_mode']}/{config['sample_prefix']}_R1_adRm.fastq.gz")
            elif config["read_mode"] == COLLAPSED:
                target_list.append(f"fastq_inputs/{config['read_mode']}/{config['sample_prefix']}_adRm.fastq.gz")
            elif config["read_mode"] == SE:
                target_list.append(f"fastq_inputs/{config['read_mode']}/{config['sample_prefix']}_adRm.fastq.gz")

        config_sample = os.path.join(args.sample_output_dir, "sample_config.yaml")

        if os.path.exists(config_sample):
            CheckExistingConfig(config_sample, config)

        with open(config_sample, "w") as fout:
            yaml.safe_dump(config, fout, default_flow_style=False)

        target_list = [os.path.join(args.sample_output_dir, target) for target in target_list]

        return self._run_snakemake("sample", args, config, target_list)

    def analyse(self):
        """
        Analyse a sample against a database
        """
        # noinspection PyTypeChecker
        parser = argparse.ArgumentParser(
            prog="haystac analyse",
            description="Analyse a sample against a database",
            formatter_class=ArgumentCustomFormatter,
            add_help=False,
        )

        required = parser.add_argument_group("Required arguments")

        required.add_argument(
            "--mode",
            help="Analysis mode for the selected sample [%(choices)s]",
            metavar="<mode>",
            choices=ANALYSIS_MODES,
            required=True,
        )

        required.add_argument(
            "--database",
            help="Path to the database output directory",
            metavar="<path>",
            dest="db_output",
            type=WritablePathType(),
            required=True,
        )

        required.add_argument(
            "--sample",
            help="Path to the sample output directory",
            metavar="<path>",
            dest="sample_output_dir",
            type=WritablePathType(),
            required=True,
        )

        required.add_argument(
            "--output",
            help="Path to the analysis output directory",
            metavar="<path>",
            dest="analysis_output_dir",
            type=WritablePathType(),
            required=True,
        )

        optional = parser.add_argument_group("Optional arguments")

        optional.add_argument(
            "--genera",
            help="List of genera to restrict the abundance calculations",
            metavar="<genus>",
            nargs="+",
            default=[],
        )

        optional.add_argument(
            "--min-prob",
            help="Minimum posterior probability to assign an aligned read to a given species (default: 0.50)",
            type=FloatRangeType(0, 100),
            metavar="<float>",
            default=argparse.SUPPRESS,
        )

        optional.add_argument(
            "--mismatch-probability",
            help="Base mismatch probability (default: 0.15)",
            type=FloatRangeType(0.01, 0.5),
            metavar="<float>",
            default=argparse.SUPPRESS,
        )

        optional.add_argument(
            "--bowtie2-threads",
            help="Number of threads bowtie2 will use to align a sample against every individual genome in the database",
            type=IntRangeType(1, MAX_CPU),
            metavar="<int>",
            dest="bowtie2_threads_aln",
            default=self.config_default["bowtie2_threads_aln"],
        )

        optional.add_argument(
            "--batch",
            help="Batch number for large`haystac analyse` workflows (e.g. --batch align_all_accessions=1/3). "
            "You will need to execute all batches before haystac is able to finish its workflow to the end.",
            metavar="<str>",
            type=BatchType(),
            default=None,
        )

        optional.add_argument(
            "--mapdamage",
            help="Perform mapdamage analysis for ancient samples (default: False)",
            type=BoolType(),
            metavar="<bool>",
            default=argparse.SUPPRESS,
        )

        optional.add_argument(
            "--aDNA",
            help="Set new flag defaults for aDNA sample analysis: "
            "mismatch-probability=0.20, min-prob=0.5, mapdamage=True",
            action="store_true",
        )

        # add the common arguments
        self._common_arguments(parser)

        # get the CLI arguments
        args = self._parse_args(parser, level=2)

        # resolve relative paths
        args.db_output = os.path.abspath(os.path.expanduser(args.db_output))
        args.sample_output_dir = os.path.abspath(os.path.expanduser(args.sample_output_dir))
        args.analysis_output_dir = os.path.abspath(os.path.expanduser(args.analysis_output_dir))

        config_fetch = os.path.join(args.db_output, "database_fetch_config.yaml")
        config_build = os.path.join(args.db_output, "database_build_config.yaml")
        config_sample = os.path.join(args.sample_output_dir, "sample_config.yaml")

        try:
            # load the database config
            database_config_file = config_fetch if os.path.exists(config_fetch) else config_build
            with open(database_config_file, "r") as fin:
                database_config = yaml.safe_load(fin)
        except FileNotFoundError:
            raise ValidationError(
                "The database has not been built correctly. Please rebuild the database using `haystac database`"
            )

        try:
            # load the sample config
            with open(config_sample, "r") as fin:
                sample_config = yaml.safe_load(fin)
        except FileNotFoundError:
            raise ValidationError(
                "The sample has not been prepared correctly. Please pre-process the samples using `haystac sample`"
            )

        # Check that the db and sample config params match
        CheckExistingConfig(database_config_file, config_sample)

        # if aDNA set new defaults
        anc_def = {"mapdamage": True, "min_prob": 0.50, "mismatch_probability": 0.20} if args.aDNA else {}

        # add all command line options to the merged config
        config = {**self.config_merged, **database_config, **sample_config, **anc_def, **vars(args)}

        # check if a db is built from an older version of haystac
        if isinstance(config["refseq_rep"], bool):
            raise ValidationError(
                "You are trying to use a database that was built with an older version of HAYSTAC. "
                "Please rebuild your database with the latest version."
            )

        target_list = list()

        if args.mode == "filter":
            if config["read_mode"] == PE:
                target_list.append(f"fastq/PE/{config['sample_prefix']}_mapq_pair.readlen")
            else:
                target_list.append(f"fastq/{config['read_mode']}/{config['sample_prefix']}_mapq.readlen")

        elif args.mode == "align":
            target_list.append(f"alignments/{config['sample_prefix']}_{config['read_mode']}_alignments.done")

        elif args.mode == "likelihoods":
            target_list.append(
                f"probabilities/{config['sample_prefix']}/{config['sample_prefix']}_likelihood_ts_tv_matrix.csv"
            )

        elif args.mode == "probabilities":
            target_list.append(
                f"probabilities/{config['sample_prefix']}/{config['sample_prefix']}_posterior_probabilities.tsv"
            )

        elif args.mode == "abundances":
            target_list.append(
                f"probabilities/{config['sample_prefix']}/{config['sample_prefix']}_posterior_abundance.tsv"
            )

        elif args.mode == "reads":
            target_list.append(f"dirichlet_reads/{config['sample_prefix']}_dirichlet_reads.done")

        # if mapdamage is true add it to the target list
        if config["mapdamage"]:
            target_list.append(f"mapdamage/{config['sample_prefix']}_mapdamage.done")

        config_analysis = os.path.join(args.analysis_output_dir, config["sample_prefix"] + "_config.yaml")

        if os.path.exists(config_analysis):
            CheckExistingConfig(config_analysis, config)

        with open(config_analysis, "w") as fout:
            yaml.safe_dump(config, fout, default_flow_style=False)

        target_list = [os.path.join(args.analysis_output_dir, target) for target in target_list]

        return self._run_snakemake("analyse", args, config, target_list)

    def _common_arguments(self, parser):
        """
        Add the common arguments shared by the `database`, `sample` and `analyse` commands
        """
        common = parser.add_argument_group("Common arguments")

        # add the help option manually so we can control where it is shown in the menu
        common.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        common.add_argument("-v", "--version", action="help", help="Print the version number and exit")

        common.add_argument(
            "--cores",
            help="Maximum number of CPU cores to use",
            metavar="<int>",
            type=IntRangeType(1, MAX_CPU),
            default=MAX_CPU if self.config_default["cores"] == "all" else self.config_default["cores"],
        )

        common.add_argument(
            "--mem",
            help="Maximum memory (MB) to use",
            type=IntRangeType(1024, MAX_MEM_MB),
            default=MAX_MEM_MB if self.config_default["mem"] == "all" else self.config_default["mem"],
            metavar="<int>",
        )

        common.add_argument(
            "--unlock",
            help="Unlock the output directory following a crash or hard restart",
            action="store_true",
        )

        common.add_argument(
            "--debug",
            help="Enable debugging mode",
            action="store_true",
        )

        common.add_argument(
            "--snakemake",
            help="Pass additional flags to the `snakemake` scheduler.",
            metavar="'<json>'",
            type=JsonType(),
        )

    def _run_snakemake(self, module, args, config, target_list):
        """
        Helper function for running the snakemake workflow
        """
        print(f"HAYSTAC v {__version__}\n")
        print(f"Date: {datetime.datetime.now()}\n")

        config["module"] = module
        config["workflow_dir"] = os.path.join(CODE_DIR, "workflow")

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

        os.makedirs(SNAKE_DIR, exist_ok=True)

        # save the run-time config file
        with open(CONFIG_RUNTIME, "w") as fout:
            yaml.safe_dump(config, fout, default_flow_style=False)

        # get any extra snakemake params
        smk_params = config.pop("snakemake") or {}

        # get any rule targets for batching and convert them into a compatible object
        exec_batch = Batch(config["batch"][0], config["batch"][1], config["batch"][2]) if config.get("batch") else None
        if config.get("batch"):
            config.pop("batch")

        success = snakemake.snakemake(
            snakefile=os.path.join(CODE_DIR, "workflow/workflow.smk"),
            config=config,
            targets=target_list,
            batch=exec_batch,
            cores=int(args.cores),
            resources={"entrez_api": self.max_entrez_requests, "mem_mb": int(args.mem)},
            force_incomplete=True,
            scheduler="greedy",
            # handle the rule-specific conda environments
            use_conda=config["use_conda"],
            conda_frontend="mamba",
            conda_prefix=os.path.join(config["cache"], "conda") if config["use_conda"] else None,
            # set all the debugging flags
            printreason=args.debug,
            printshellcmds=args.debug,
            show_failed_logs=args.debug,
            verbose=args.debug,
            keep_incomplete=args.debug,
            restart_times=0 if args.debug else RESTART_TIMES,
            keepgoing=(not args.debug),
            unlock=args.unlock,
            # pass on any CLI arguments from the --snakemake flag
            **smk_params,
        )

        # tidy up all the snakemake metadata
        if success and (len(os.listdir(os.path.join(SNAKE_DIR, "locks"))) == 0):
            shutil.rmtree(SNAKE_DIR)

        # translate "success" into shell exit code of 0
        return 0 if success else 1


if __name__ == "__main__":
    Haystac()
