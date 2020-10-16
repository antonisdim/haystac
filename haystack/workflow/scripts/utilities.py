#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import argparse
import os
import pandas as pd
import re
import sys
import time
import urllib.error
from Bio import Entrez

sys.path.append(os.getcwd())
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from scripts.entrez_utils import ENTREZ_DB_ASSEMBLY

TOO_MANY_REQUESTS_WAIT = 20
MAX_RETRY_ATTEMPTS = 5


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


class EmailType(object):
    """
    Supports checking email against different patterns. The current available patterns is:
    RFC5322 (http://www.ietf.org/rfc/rfc5322.txt)

    See https://gist.github.com/asfaltboy/79a02a2b9871501af5f00c95daaeb6e7
    """

    patterns = {
        "RFC5322": re.compile(r"^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$"),
    }

    def __init__(self, pattern):
        if pattern not in self.patterns:
            raise KeyError(
                "{} is not a supported email pattern, choose from:" " {}".format(pattern, ",".join(self.patterns))
            )
        self._rules = pattern
        self._pattern = self.patterns[pattern]

    def __call__(self, value):
        if not self._pattern.match(value):
            raise argparse.ArgumentTypeError(f"'{value}' is not a valid email - does not match {self._rules} rules")
        return value


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
                f"'{value}' is not a valid {self.type.__name__} in the range " f"({self.lower}, {self.upper})"
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


def get_total_paths(
    wildcards, checkpoints, entrez_query, with_refseq_rep, sequences, accessions, specific_genera,
):
    """
    Get all the individual fasta file paths for the taxa in our database.
    """

    sequences_df = pd.DataFrame()

    if entrez_query:
        pick_sequences = checkpoints.entrez_pick_sequences.get()
        sequences_df = pd.read_csv(pick_sequences.output[0], sep="\t")

        if len(sequences_df) == 0:
            raise RuntimeError("The entrez pick sequences file is empty.")

    if with_refseq_rep:
        refseq_rep_prok = checkpoints.entrez_refseq_accessions.get()

        refseq_genomes = pd.read_csv(refseq_rep_prok.output.refseq_genomes, sep="\t")
        genbank_genomes = pd.read_csv(refseq_rep_prok.output.genbank_genomes, sep="\t")
        assemblies = pd.read_csv(refseq_rep_prok.output.assemblies, sep="\t")
        refseq_plasmids = pd.read_csv(refseq_rep_prok.output.refseq_plasmids, sep="\t")
        genbank_plasmids = pd.read_csv(refseq_rep_prok.output.genbank_plasmids, sep="\t")

        invalid_assemblies = checkpoints.entrez_invalid_assemblies.get()
        invalid_assembly_sequences = pd.read_csv(invalid_assemblies.output[0], sep="\t")

        assemblies = assemblies[
            ~assemblies["GBSeq_accession-version"].isin(invalid_assembly_sequences["GBSeq_accession-version"])
        ]

        sources = [
            refseq_genomes,
            genbank_genomes,
            assemblies,
            refseq_plasmids,
            genbank_plasmids,
        ]

        if entrez_query:
            sources.append(sequences_df)

        sequences_df = pd.concat(sources)

    if sequences:
        custom_fasta_paths = pd.read_csv(
            sequences, sep="\t", header=None, names=["species", "GBSeq_accession-version", "path"],
        )

        custom_seqs = custom_fasta_paths[["species", "GBSeq_accession-version"]]
        custom_seqs["GBSeq_accession-version"] = "custom_seq-" + custom_seqs["GBSeq_accession-version"].astype(str)

        sequences_df = sequences_df.append(custom_seqs)

    if accessions:
        custom_accessions = pd.read_csv(
            accessions, sep="\t", header=None, names=["species", "GBSeq_accession-version"],
        )

        sequences_df = sequences_df.append(custom_accessions)

    if specific_genera:
        sequences_df = sequences_df[sequences_df["species"].str.contains("|".join(specific_genera))]

    return sequences_df


def normalise_name(taxon):
    """remove unnecessary characters from a taxon name string."""
    return taxon.replace(" ", "_").replace("[", "").replace("]", "").replace("/", "_")


def check_unique_taxa_in_custom_input(accessions, sequences):
    """Checks that custom input files have only one entry per taxon"""

    if accessions != "" and sequences != "":
        custom_fasta_paths = pd.read_csv(
            config["sequences"], sep="\t", header=None, names=["species", "accession", "path"],
        )
        custom_accessions = pd.read_csv(config["accessions"], sep="\t", header=None, names=["species", "accession"],)

        taxon_acc = custom_accessions["species"].tolist()
        taxon_seq = custom_fasta_paths["species"].tolist()

        if bool(set(taxon_acc) & set(taxon_seq)):
            raise RuntimeError(
                "You have provided the same taxon both in your custom sequences file and your "
                "custom accessions file. Please pick and keep ONLY one entry from both of these files. "
                "You can only have 1 sequence per chosen taxon in your database."
            )


def get_accession_ftp_path(accession, config, attempt=1):
    """Get a valid NCBI ftp path from an accession."""

    Entrez.email = config["email"]
    try:
        handle = Entrez.esearch(db=ENTREZ_DB_ASSEMBLY, term=accession + ' AND "latest refseq"[filter]')
        # or handle = Entrez.esearch(db=ENTREZ_DB_ASSEMBLY,
        # term=accession + ' AND ((latest[filter] OR "latest refseq"[filter])')
        assembly_record = Entrez.read(handle)
        esummary_handle = Entrez.esummary(db=ENTREZ_DB_ASSEMBLY, id=assembly_record["IdList"], report="full")
        try:
            esummary_record = Entrez.read(esummary_handle, validate=False)
        except RuntimeError:
            return ""
        refseq_ftp = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
        genbank_ftp = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_GenBank"]

        if refseq_ftp != "":
            return refseq_ftp
        else:
            return genbank_ftp

    except urllib.error.HTTPError as e:
        if e.code == 429:

            attempt += 1

            if attempt > MAX_RETRY_ATTEMPTS:
                print("Exceeded maximum attempts {}...".format(attempt), file=sys.stderr)
                return None
            else:
                time.sleep(TOO_MANY_REQUESTS_WAIT)
                get_accession_ftp_path(accession, config, attempt)

        else:
            raise RuntimeError("There was a urllib.error.HTTPError with code {}".format(e))

    except IndexError:
        time.sleep(TOO_MANY_REQUESTS_WAIT)
        get_accession_ftp_path(accession, config, attempt)
