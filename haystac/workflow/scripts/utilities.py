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

WARNING_DB = 0
WARNING_USER = 0

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
            self.data = pd.read_table(
                value,
                sep="\t",
                header=None,
                index_col=False,
            )
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

        run_code = etree.find(".//RUN").attrib["accession"]

        if len(id_list) > 1 or value != run_code:
            raise argparse.ArgumentTypeError(
                f"The SRA accession you have provided {value} does not refer to a sequencing run. "
                f"Please visit https://www.ncbi.nlm.nih.gov/sra/ and chose a valid "
                f"sequencing run accession for the SRA accession {value}."
            )

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
                        "aDNA",
                        "seed",
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


class BatchType(object):
    """
    Is this a valid smk batch string
    """

    def __call__(self, value):

        try:
            rulename, batch, batches = (
                value.split("=")[0],
                int(value.split("=")[1].split("/")[0]),
                int(value.split("=")[1].split("/")[1]),
            )
            return rulename, batch, batches
        except IndexError as error:
            raise argparse.ArgumentTypeError(f"'{value}' is not a valid snakemake batch string\n {error}")


def get_total_paths(
    checkpoints,
    config,
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

        sources = []

        # refseq rep prok
        if config["refseq_rep"] == "prokaryote_rep":
            refseq_rep_prok = checkpoints.entrez_refseq_rep_prok_accessions.get()
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

        # refseq viruses
        elif config["refseq_rep"] == "viruses":
            refseq_viruses = checkpoints.entrez_refseq_viruses_accessions.get()
            refseq_viral_genomes = pd.read_csv(refseq_viruses.output.refseq_viruses, sep="\t")

            sources = [refseq_viral_genomes]

        # refseq eukaryotes
        elif config["refseq_rep"] == "eukaryotes":
            refseq_eukaryotes = checkpoints.entrez_refseq_eukaryotes_accessions.get()
            refseq_euk_genomes = pd.read_csv(refseq_eukaryotes.output.refseq_euk, sep="\t")

            sources = [refseq_euk_genomes]

        if config["query"]:
            sources.append(sequences_df)

        sequences_df = pd.concat(sources)

    if config["sequences"] or config["accessions"]:
        check_unique_taxa_in_custom_inputs(config["accessions"], config["sequences"])

    if config["sequences"]:
        custom_fasta_paths = pd.read_csv(
            config["sequences"],
            sep="\t",
            header=None,
            names=["species", "AccessionVersion", "path"],
        )
        custom_fasta_paths = check_unique_taxa_accs(custom_fasta_paths, config, config["sequences"], "user_file")

        custom_seqs = custom_fasta_paths[["species", "AccessionVersion"]].copy()
        custom_seqs["AccessionVersion"] = "custom_seq-" + custom_seqs["AccessionVersion"].astype(str)

        sequences_df = sequences_df.append(custom_seqs)

    if config["accessions"]:
        custom_accessions = pd.read_csv(
            config["accessions"],
            sep="\t",
            header=None,
            names=["species", "AccessionVersion"],
        )
        custom_accessions = check_unique_taxa_accs(custom_accessions, config, config["accessions"], "user_file")

        sequences_df = sequences_df.append(custom_accessions)

    if config["genera"]:
        sequences_df = sequences_df[sequences_df["species"].str.contains("|".join(config["genera"]))]

    if config["exclude_accessions"]:
        sequences_df = sequences_df[~sequences_df["AccessionVersion"].isin(config["exclude_accessions"])]

    # check that db accessions are unique
    sequences_df = check_unique_taxa_accs(sequences_df, config, "", "db")

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


def check_unique_taxa_in_custom_inputs(accessions, sequences):
    """Checks that custom input files have only one entry per taxon"""

    if accessions != "" and sequences != "":
        custom_fasta_paths = pd.read_csv(sequences, sep="\t", header=None, names=["species", "accession", "path"])
        custom_accessions = pd.read_csv(accessions, sep="\t", header=None, names=["species", "accession"])

        # check if any taxa in common
        taxon_acc = custom_accessions["species"].tolist()
        taxon_seq = custom_fasta_paths["species"].tolist()

        if bool(set(taxon_acc) & set(taxon_seq)):
            print_error(
                "You have provided the same taxon both in your custom sequences "
                "file and your custom accessions file. Please pick and keep ONLY "
                "one entry from both of these files. You can only have 1 sequence "
                "per chosen taxon in your database."
            )

        # check if any accessions in common
        accession_acc = custom_accessions["accession"].tolist()
        accession_seq = custom_fasta_paths["accession"].tolist()

        if bool(set(accession_acc) & set(accession_seq)):
            print_error(
                "You have provided the same accession both in your custom sequences "
                "file and your custom accessions file. Please pick and keep ONLY "
                "one entry from both of these files, or change the accession entry "
                "appropriately in your custom sequences file. You can only have 1 accession name "
                "per chosen taxon in your database."
            )


def check_unique_taxa_accs(df, config, user_input, to_check):
    """Checks that there are only unique inputs for taxa and accessions"""

    # if we are checking the user files
    if to_check == "user_file":

        # if duplicate accession in user file raise error
        if df["AccessionVersion"].duplicated().any():
            dup_acc = [i for i in df[df["AccessionVersion"].duplicated()]["AccessionVersion"].to_list()]
            message = (
                f"{user_input} contains multiple taxa for {', '.join(dup_acc)}. "
                f"Please remove/fix all duplicates. Picking automatically a taxon/accession pair in "
                f"this case is not possible."
            )
            print_error(message)

        # if duplicate species in user file either raise error, or --resolve-accessions
        elif df["species"].duplicated().any():
            dup_taxa = [i for i in df[df["species"].duplicated()]["species"].to_list()]
            message = f"{user_input} contains multiple sequences for {', '.join(dup_taxa)}. "

            if not config["resolve_accessions"]:
                message += (
                    "Either remove all duplicates, or set the `--resolve-accessions` flag to automatically choose one. "
                    "It is the first accession that will be chosen."
                )
                print_error(message)
            else:
                # global WARNING_USER
                for idx, val in df[df["species"].duplicated(keep="first")].iterrows():
                    message += f"Accession {val['AccessionVersion']} for {val['species']} was omitted."
                    # if WARNING_USER == 0:
                    print_warning(message)
                # WARNING_USER += 1
                df = df[~df["species"].duplicated(keep="first")]
                return df

        # if all good return df as is
        else:
            return df

    # if we are checking the database in total
    elif to_check == "db":

        # if duplicate accessions in db either raise error, or --resolve-accessions
        if df["AccessionVersion"].duplicated().any():
            dup_acc = [i for i in df[df["AccessionVersion"].duplicated()]["AccessionVersion"].to_list()]
            dup_tax = [i for i in df[df["AccessionVersion"].duplicated(keep=False)]["species"].to_list()]
            message = (
                f"Accession {', '.join(dup_acc)} appears multiple times in the database "
                f"with different taxa names ({', '.join(dup_tax)}). "
            )

            if not config["resolve_accessions"]:
                message += (
                    f"Please remove/fix all duplicate accessions if possible. "
                    f"If multiple taxa have the same accession, "
                    f"that is possibly due to a recent change in NCBI's taxonomy, and it is strongly "
                    f"advised you check the latest information for these accessions. "
                    f"Either specify unique pairs of taxa and accessions using the `--accessions-file` or "
                    f"`--sequences-file` flags, or set the `--resolve-accessions` flag to automatically "
                    f"choose the first one. "
                )
                print_error(message)
            else:
                # global WARNING_DB
                for idx, val in df[df["AccessionVersion"].duplicated(keep="first")].iterrows():
                    message += (
                        f"{val['species']} has been excluded. It is strongly advised to "
                        f"check the latest taxonomy info on NCBI."
                    )
                    # if WARNING_DB == 0:
                    print_warning(message)
                # WARNING_DB += 1
                df = df[~df["AccessionVersion"].duplicated(keep="first")]
                return df

        # if all good return df as id
        else:
            return df


def get_final_db_paths(checkpoints):
    """Get all the taxon/acc pairs for the taxa in our database."""

    db_sequences = checkpoints.entrez_db_list.get()
    sequences_df = pd.read_csv(db_sequences.output[0], sep="\t", names=["species", "AccessionVersion"])

    assert len(sequences_df) > 0, (
        f"The db file containing the taxon/accession pairs is empty {db_sequences.output[0]}. "
        f"Please rebuild the database."
    )

    inputs = []

    for key, seq in sequences_df.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["AccessionVersion"],
        )
        inputs.append((orgname, accession))

    return inputs


def chunker(seq, size):
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


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


def get_smk_config():
    """Function to read the smk config and return a dictionary"""

    try:
        with open(".snakemake/config.yaml") as fin:
            config = yaml.safe_load(fin)
    except FileNotFoundError:
        config = {}

    return config
