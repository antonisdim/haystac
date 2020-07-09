#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import urllib.error
import time

from Bio import Entrez
from Bio import bgzf

sys.path.append(os.getcwd())
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from scripts.entrez_utils import (
    guts_of_entrez,
    ENTREZ_DB_NUCCORE,
    ENTREZ_RETMODE_TEXT,
    ENTREZ_RETTYPE_FASTA,
    ENTREZ_DB_ASSEMBLY,
    ENTREZ_RETMODE_XML,
    ENTREZ_RETTYPE_GB,
)

TOO_MANY_REQUESTS_WAIT = 5
MAX_RETRY_ATTEMPTS = 3

# TODO refactor this code to make it simpler and cleaner
def entrez_download_sequence(accession, config, output_file, attempt=1, assembly=False):
    """
    Fetch the reference genome from NCBI.
    """
    print("The sequences for the database are being downloaded ...", file=sys.stderr)

    Entrez.email = config["entrez_email"]

    if assembly:
        handle = Entrez.esearch(
            db=ENTREZ_DB_ASSEMBLY, term=accession + ' AND "latest refseq"[filter]'
        )
        # or handle = Entrez.esearch(db=ENTREZ_DB_ASSEMBLY,
        # term=accession + ' AND ((latest[filter] OR "latest refseq"[filter])')
        assembly_record = Entrez.read(handle)
        esummary_handle = Entrez.esummary(
            db=ENTREZ_DB_ASSEMBLY, id=assembly_record["IdList"], report="full"
        )
        esummary_record = Entrez.read(esummary_handle)
        accession_id = esummary_record["DocumentSummarySet"]["DocumentSummary"][0][
            "AssemblyAccession"
        ]
        nuccore_id = Entrez.read(
            Entrez.esearch(db=ENTREZ_DB_NUCCORE, term=accession_id)
        )["IdList"]
    else:
        nuccore_id = [accession]

    try:
        records = guts_of_entrez(
            ENTREZ_DB_NUCCORE,
            ENTREZ_RETMODE_TEXT,
            ENTREZ_RETTYPE_FASTA,
            nuccore_id,
            batch_size=1,
        )

        if all(
            elem == "\n"
            for elem in [
                fa
                for fa in guts_of_entrez(
                    ENTREZ_DB_NUCCORE,
                    ENTREZ_RETMODE_TEXT,
                    ENTREZ_RETTYPE_FASTA,
                    nuccore_id,
                    batch_size=1,
                )
            ]
        ):
            print(
                "The accession {} is a master record for an WGS project. "
                "That means that this record is empty in nuccore. "
                "Going to the assembly database to fetch the assembly for this taxon.".format(
                    accession
                ),
                file=sys.stderr,
            )
            xml_records = guts_of_entrez(
                ENTREZ_DB_NUCCORE,
                ENTREZ_RETMODE_XML,
                ENTREZ_RETTYPE_GB,
                nuccore_id,
                batch_size=1,
            )
            xml_list = [xml for xml in xml_records]
            assembly_id = xml_list[0]["GBSeq_xrefs"][2]["GBXref_id"]
            new_nuccore_id = Entrez.read(
                Entrez.esearch(db=ENTREZ_DB_NUCCORE, term=assembly_id)
            )["IdList"]
            records = guts_of_entrez(
                ENTREZ_DB_NUCCORE,
                ENTREZ_RETMODE_TEXT,
                ENTREZ_RETTYPE_FASTA,
                new_nuccore_id,
                batch_size=1,
            )

        with bgzf.open(output_file, "wt") as fout:
            for fasta in records:
                fout.write(fasta)

    except urllib.error.HTTPError as e:
        if e.code == 429:

            attempt += 1

            if attempt > MAX_RETRY_ATTEMPTS:
                print(
                    "Exceeded maximum attempts {}...".format(attempt), file=sys.stderr
                )
                return None
            else:
                time.sleep(TOO_MANY_REQUESTS_WAIT)
                entrez_download_sequence(accession, config, output_file, attempt)

        else:
            raise RuntimeError(
                "There was a urllib.error.HTTPError with code {}".format(e)
            )


if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    entrez_download_sequence(
        accession=snakemake.wildcards.accession,
        config=snakemake.config,
        output_file=snakemake.output[0],
        assembly=snakemake.params[0],
    )
