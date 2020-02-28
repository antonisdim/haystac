#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from pprint import pprint

from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import ENTREZ_DB_NUCCORE, guts_of_entrez, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_GB

Entrez.email = 'noreply@example.com'
entrez_query = 'NZ_CP019458.1[Accession] OR NZ_CP021748.1[Accession] OR NZ_LT670818.1[Accession] ' \
               'OR NZ_CP022415.1[Accession] OR NZ_FOYO01000001.1[Accession] OR NZ_CP020121.1[Accession]'

retmax = 2
counter = 0


def gen_dict_extract(key, var):
    """Find keys in nested dictionaries"""
    if hasattr(var, 'items'):
        for k, v in var.items():
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result
    elif isinstance(var, list):
        for d in var:
            for result in gen_dict_extract(key, d):
                yield result


while True:
    handle = Entrez.esearch(db=ENTREZ_DB_NUCCORE, term=entrez_query, retmax=retmax, idtype="acc", usehistory='Y',
                            retstart=retmax * counter, rettype='gb', retmode='xml')
    handle = Entrez.read(handle)

    if not handle['IdList']:
        # stop iterating when we get an empty resultset
        break

    # TODO set the following to save downloading the full sequence each time
    records = list(guts_of_entrez(ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_GB, handle['IdList'], retmax))

    for record in records:
        # extract the taxon id from the overly complicated nesting in the GBSeq XML format
        record['TSeq_taxid'] = [val.replace('taxon:', '')
                                for val in gen_dict_extract('GBQualifier_value', record) if 'taxon:' in val].pop()

    pprint(records)

    counter += 1

    break
