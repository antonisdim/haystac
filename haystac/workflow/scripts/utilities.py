#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import argparse
import hashlib
import os
import re
import sys

import pandas as pd
import yaml

REGEX_WHITELIST = r"[\w.-]+"
REGEX_BLACKLIST = r"[^\w.-]+"

PE = "PE"
COLLAPSED = "COLLAPSED"
SE = "SE"

WARNING = "\x1b[33m"
FAIL = "\x1b[31m"
END = "\033[0m"

is_tty = sys.stdout.isatty()


class ValidationError(Exception):
    pass


class ArgumentCustomFormatter(argparse.HelpFormatter):
    """
    Custom formatter for argparse
    """

    def _get_help_string(self, action):
        message = action.help
        if "%(default)" not in action.help:
            if action.default is not argparse.SUPPRESS and action.default is not None:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    message += " (default: %(default)s)"
        return message


class FileType(argparse.FileType):
    """
    Override argparse.FileType to return the filename, rather than an open file handle.
    """

    def __call__(self, string):
        return super().__call__(string).name


class WritablePathType(object):
    """
    Is this a writable path.
    """

    def __call__(self, value):
        from pathlib import Path

        try:
            path = Path(value).expanduser()
            path.mkdir(parents=True, exist_ok=True)
            return value
        except Exception:
            raise argparse.ArgumentTypeError(f"'{value}' is not a valid writable path")


class PositiveIntType(object):
    """
    Is this a positive integer
    """

    def __call__(self, value):
        try:
            if not int(value) > 0:
                raise ValueError()
        except ValueError:
            raise argparse.ArgumentTypeError(f"'{value}' is not a valid positive integer")

        return int(value)


class RangeType(object):
    """
    Is this a valid instance of `_type` and within the range [lower, upper]
    """

    def __init__(self, _type, lower, upper):
        self.type = _type
        self.lower = lower
        self.upper = upper

    def __call__(self, value):
        try:
            if not (self.lower <= self.type(value) <= self.upper):
                raise ValueError()
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"'{value}' is not a valid {self.type.__name__} in the range ({self.lower}, {self.upper})"
            )

        return self.type(value)


class FloatRangeType(RangeType):
    """
    Is this a float() within the given range
    """

    def __init__(self, lower, upper):
        super().__init__(float, lower, upper)


class IntRangeType(RangeType):
    """
    Is this an int() within the given range
    """

    def __init__(self, lower, upper):
        super().__init__(int, lower, upper)


class BoolType(object):
    """
    Is this a valid boolean
    """

    def __call__(self, value):
        if isinstance(value, bool):
            return value
        if value.lower() in ("yes", "true", "t", "y", "1"):
            return True
        elif value.lower() in ("no", "false", "f", "n", "0"):
            return False
        else:
            raise argparse.ArgumentTypeError(f"'{value}' is not a valid boolean")


class JsonType(object):
    """
    Is this a valid JSON string
    """

    def __call__(self, value):
        import json

        try:
            return json.loads(value)
        except json.decoder.JSONDecodeError as error:
            raise argparse.ArgumentTypeError(f"'{value}' is not a valid JSON string\n {error}")


class SpreadsheetFileType(object):
    """
    Is it a valid user input file
    """

    cols = None
    data = None

    def __call__(self, value):

        if not os.path.exists(value):
            raise argparse.ArgumentTypeError(f"'{value}' does not exit")

        if os.stat(value).st_size == 0:
            raise argparse.ArgumentTypeError(f"'{value}' is empty")

        try:
            self.data = pd.read_table(value, sep="\t", header=None, index_col=False,)
        except Exception:
            raise argparse.ArgumentTypeError(f"'{value}' unknown error parsing file")

        if len(self.data.columns) != len(self.cols):
            raise argparse.ArgumentTypeError(
                f"'{value}' must have {len(self.cols)} columns and be tab delimited (cols={len(self.data.columns)})"
            )

        # find the row number of any empty cells
        bad_rows = ", ".join([str(i + 1) for i in self.data.index[self.data.isnull().any(axis=1)].tolist()])

        if bad_rows:
            raise argparse.ArgumentTypeError(f"'{value}' contains missing data in line(s): {bad_rows}")

        return value


class AccessionFileType(SpreadsheetFileType):
    """
    Is this a valid accession input file.
    """

    cols = ["species", "accession"]

    def __call__(self, value):
        super().__call__(value)

        # check all accessions pass the regex pattern
        idx = self.cols.index("accession")
        species = self.cols.index("species")

        bad_list_acc = "\n".join(
            [
                f"line {i + 1}: '{acc}'"
                for i, acc in enumerate(self.data[idx].tolist())
                if re.match(REGEX_BLACKLIST, acc) is not None
            ]
        )

        bad_list_species = "\n".join(
            [
                f"line {i + 1}: '{tax}'"
                for i, tax in enumerate(self.data[species].tolist())
                if re.match(REGEX_BLACKLIST, tax) is not None
            ]
        )

        if bad_list_acc or bad_list_species:
            bad_list = bad_list_acc + "\n" + bad_list_species
            raise argparse.ArgumentTypeError(
                f"'{value}' these accession codes or taxon names contain invalid characters:\n{bad_list}"
            )

        return value


class SequenceFileType(AccessionFileType):
    """
    Is this a valid sequence input file.
    """

    cols = ["species", "accession", "path"]

    def __call__(self, value):
        super().__call__(value)

        # find any files that don't exist or are empty
        idx = self.cols.index("path")

        bad_files = "\n".join(
            [
                f"line {i + 1}: '{file}'"
                for i, file in enumerate(self.data[idx].tolist())
                if not os.path.exists(file) or os.stat(file).st_size == 0
            ]
        )

        if bad_files:
            raise argparse.ArgumentTypeError(f"'{value}' these sequence files do not exist or are empty:\n{bad_files}")

        return value


class SraAccessionType(object):
    """
    Is this a valid SRA accession
    """

    def __call__(self, value):
        # import these locally to avoid cyclic import issues
        from haystac.workflow.scripts.entrez_utils import entrez_esearch, entrez_efetch

        try:
            # query the SRA to see if this a valid accession
            _, _, id_list = entrez_esearch("sra", f"{value}[Accession]")
            etree = entrez_efetch("sra", id_list)
        except Exception:
            raise argparse.ArgumentTypeError(f"Invalid SRA accession '{value}'")

        try:
            # now get the library layout
            layout = etree.find(".//LIBRARY_LAYOUT/*").tag.lower()
        except Exception:
            raise argparse.ArgumentTypeError(f"Unable to resolve the library layout for SRA accession '{value}'")

        return value, layout


class NuccoreQueryType(object):
    """
    Is this a valid nuccore query
    """

    def __call__(self, value):
        # import these locally to avoid cyclic import issues
        from haystac.workflow.scripts.entrez_utils import entrez_esearch

        # check if the user has given us a file instead of a string
        if os.path.isfile(value):
            query = open(value).read().strip()
            if not query:
                raise argparse.ArgumentTypeError(f"The query file '{value}' is empty.")
        else:
            query = value

        try:
            # query nuccore to see if this a valid query
            _, _, id_list = entrez_esearch("nuccore", f"{query}")
        except Exception:
            raise argparse.ArgumentTypeError(f"Invalid NCBI query '{query}'")

        # if the query returns no result set raise error
        if len(id_list) == 0:
            raise argparse.ArgumentTypeError(f"No results in NCBI nucleotide for query '{query}'")

        return value


class CheckExistingConfig(object):
    """
    Checks the details of an existing yaml file against cli params or another yaml file
    """

    def __init__(self, filename, params):

        # check if second argument is a dict or a yaml file
        if isinstance(params, dict):
            params_config = params
        elif os.path.isfile(params):
            with open(params, "r") as fin_params:
                params_config = yaml.safe_load(fin_params)

        # check if a config already exists
        if not os.path.isfile(filename):
            pass

        else:
            # open the config file
            with open(filename, "r") as fin:
                existing_config = yaml.safe_load(fin)

                if not isinstance(params, dict):
                    important_args = ["cache"]
                else:
                    important_args = [
                        "cache",
                        "api_key",
                        "mismatch_probability",
                        "bowtie2_scaling",
                        "query",
                        "query_file",
                        "accessions_file",
                        "sequences_file",
                        "refseq_rep",
                        # "force_accessions",
                        # "exclude_accessions",
                        # "resolve_accessions",
                        "rank",
                        "mtDNA",
                        "seed",
                        "sample_prefix",
                        "fastq",
                        "fastq_r1",
                        "fastq_r2",
                        "sra",
                        "collapse",
                        "trim_adapters",
                        "sample",
                        "min_prob",
                        "query_file_md5",
                        "accessions_md5",
                        "sequences_md5",
                    ]

                for arg in important_args:
                    # check if all the important params match
                    if arg in existing_config.keys() and arg in params_config.keys():
                        if existing_config[arg] != params_config[arg]:
                            print_error(
                                f"You are trying to set a value for parameter {arg} on top of an already existing one "
                                f"(old: {existing_config[arg]}, new: {params_config[arg]}). "
                                f"Please either revert to the original parameter you used or create a "
                                f"new output directory."
                            )


class FastqFile(object):
    """
    Is it a valid user input fastq file
    """

    def __call__(self, value):

        if not os.path.exists(value):
            raise argparse.ArgumentTypeError(f"'{value}' does not exit")

        if os.stat(value).st_size == 0:
            raise argparse.ArgumentTypeError(f"'{value}' is empty")

        if ".gz" not in value:
            with open(value, "r") as fin:
                first_line = fin.readline()
        else:
            import gzip

            with gzip.open(value, "rt") as fin:
                first_line = fin.readline()

        if first_line[0] != "@":
            raise argparse.ArgumentTypeError(f"'{value}' is not a valid fastq file.")

        return value


def get_total_paths(
    checkpoints, config,
):
    """
    Get all the individual fasta file paths for the taxa in our database.
    """

    sequences_df = pd.DataFrame()

    if config["query"]:
        pick_sequences = checkpoints.entrez_pick_sequences.get()
        sequences_df = pd.read_csv(pick_sequences.output[0], sep="\t")

        assert len(sequences_df) > 0, f"The entrez pick sequences file is empty {pick_sequences.output[0]}"

    if config["refseq_rep"]:
        refseq_rep_prok = checkpoints.entrez_refseq_accessions.get()
        refseq_genomes = pd.read_csv(refseq_rep_prok.output.refseq_genomes, sep="\t")
        genbank_genomes = pd.read_csv(refseq_rep_prok.output.genbank_genomes, sep="\t")
        assemblies = pd.read_csv(refseq_rep_prok.output.assemblies, sep="\t")
        refseq_plasmids = pd.read_csv(refseq_rep_prok.output.refseq_plasmids, sep="\t")
        genbank_plasmids = pd.read_csv(refseq_rep_prok.output.genbank_plasmids, sep="\t")

        if not config["force_accessions"]:
            invalid_assemblies = checkpoints.entrez_invalid_assemblies.get()
            invalid_assembly_sequences = pd.read_csv(invalid_assemblies.output[0], sep="\t")

            assemblies = assemblies[
                ~assemblies["AccessionVersion"].isin(invalid_assembly_sequences["AccessionVersion"])
            ]

        sources = [
            refseq_genomes,
            genbank_genomes,
            assemblies,
            refseq_plasmids,
            genbank_plasmids,
        ]

        if config["query"]:
            sources.append(sequences_df)

        sequences_df = pd.concat(sources)

    if config["sequences"] or config["accessions"]:
        check_unique_taxa_in_custom_inputs(config["accessions"], config["sequences"])

    if config["sequences"]:
        custom_fasta_paths = pd.read_csv(
            config["sequences"], sep="\t", header=None, names=["species", "AccessionVersion", "path"],
        )
        check_unique_taxa_in_user_file(custom_fasta_paths, config)

        custom_seqs = custom_fasta_paths[["species", "AccessionVersion"]].copy()
        custom_seqs["AccessionVersion"] = "custom_seq-" + custom_seqs["AccessionVersion"].astype(str)

        sequences_df = sequences_df.append(custom_seqs)

    if config["accessions"]:
        custom_accessions = pd.read_csv(
            config["accessions"], sep="\t", header=None, names=["species", "AccessionVersion"],
        )
        check_unique_taxa_in_user_file(custom_accessions, config)

        sequences_df = sequences_df.append(custom_accessions)

    if config["genera"]:
        sequences_df = sequences_df[sequences_df["species"].str.contains("|".join(config["genera"]))]

    if config["exclude_accessions"]:
        sequences_df = sequences_df[~sequences_df["AccessionVersion"].isin(config["exclude_accessions"])]

    inputs = []

    for key, seq in sequences_df.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["AccessionVersion"],
        )
        inputs.append((orgname, accession))

    return inputs


def normalise_name(taxon):
    """remove unnecessary characters from a taxon name string."""

    return re.sub(REGEX_BLACKLIST, "_", taxon)


# return taxon.replace(" ", "_").replace("[", "").replace("]", "").replace("/", "_")


def check_unique_taxa_in_custom_inputs(accessions, sequences):
    """Checks that custom input files have only one entry per taxon"""

    if accessions != "" and sequences != "":
        custom_fasta_paths = pd.read_csv(sequences, sep="\t", header=None, names=["species", "accession", "path"])
        custom_accessions = pd.read_csv(accessions, sep="\t", header=None, names=["species", "accession"])

        taxon_acc = custom_accessions["species"].tolist()
        taxon_seq = custom_fasta_paths["species"].tolist()

        if bool(set(taxon_acc) & set(taxon_seq)):
            print_error(
                "You have provided the same taxon both in your custom sequences "
                "file and your custom accessions file. Please pick and keep ONLY "
                "one entry from both of these files. You can only have 1 sequence "
                "per chosen taxon in your database."
            )


def chunker(seq, size):
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


def check_unique_taxa_in_user_file(df, config):
    """Checks that there are only unique inputs for taxa"""

    if df["species"].duplicated().any():
        dup_taxa = [i for i in df[df["species"].duplicated()]["species"].to_list()]
        message = f"{config['sequences']} contains multiple sequences for {', '.join(dup_taxa)}. "

        if not config["resolve_accessions"]:
            message += (
                "Either remove all duplicates, or set the `--resolve-accessions` flag to automatically choose one. "
                "It is the first accession that will be chosen."
            )
            print_error(message)
        else:
            df = df[~df["species"].duplicated(keep="first")]
            for idx, val in df[df["species"].duplicated()].iterrows():
                message += f"Accession {val['accession']} for {val['species']} was omitted."
                print_warning(message)
            return df


def md5(filename):
    hash_md5 = hashlib.md5()

    # open file and get the checksum
    with open(filename, "rb") as f:
        # read it in chunks in case the file is big
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)

    return hash_md5.hexdigest()


def print_error(message):
    """Function to print errors and exit"""
    message = f"haystac: error: {message}"
    print(f"{FAIL}{message}{END}" if is_tty else message, file=sys.stderr)
    exit(1)


def print_warning(message):
    """Function to print warnings"""
    message = f"WARNING: {message}"
    print(f"{WARNING}{message}{END}" if is_tty else message, file=sys.stderr)
