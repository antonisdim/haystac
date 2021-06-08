#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from collections import defaultdict

import pysam


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):

        if not read.is_proper_pair:
            continue

        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def mutation_type(alleles):
    """
    Is this mutation a transition (A <-> G and C <-> T) or a transversion (everything else).
    """
    alleles = set([a.upper() for a in alleles])
    assert len(alleles) == 2 and alleles.issubset({"A", "C", "G", "T", "N"})
    return "ts" if alleles in [{"A", "G"}, {"C", "T"}] else "tv"


def count_ts_tv_init(bam_file, output_file, taxon, pairs=False):
    bam = pysam.AlignmentFile(bam_file, "rb")

    with open(output_file, "w") as fout:

        for read in bam.fetch():
            ts, tv = 0, 0

            if not read.is_proper_pair:

                if pairs:
                    continue

                for base_call, base_ref in zip(read.seq, read.get_reference_sequence()):
                    if base_call.upper() == "N" or base_ref.upper() == "N":
                        continue

                    if base_call != base_ref:
                        if mutation_type([base_call, base_ref]) == "ts":
                            ts += 1
                        else:
                            tv += 1

                print(taxon, read.query_name, ts, tv, file=fout, sep=",")

        if pairs:

            for read1, read2 in read_pair_generator(bam):
                ts, tv = 0, 0

                for base_call1, base_ref1, base_call2, base_ref2 in zip(
                    read1.seq,
                    read1.get_reference_sequence(),
                    read2.seq,
                    read2.get_reference_sequence(),
                ):

                    if (
                        base_call1.upper() == "N"
                        or base_ref1.upper() == "N"
                        or base_call2.upper() == "N"
                        or base_ref2.upper() == "N"
                    ):
                        continue

                    if base_call1 != base_ref1:
                        if mutation_type([base_call1, base_ref1]) == "ts":
                            ts += 1
                        else:
                            tv += 1

                    if base_call2 != base_ref2:
                        if mutation_type([base_call2, base_ref2]) == "ts":
                            ts += 1
                        else:
                            tv += 1

                print(taxon, read1.query_name, ts, tv, file=fout, sep=",")


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    count_ts_tv_init(
        bam_file=snakemake.input[0],
        output_file=snakemake.output[0],
        taxon=snakemake.wildcards.orgname,
        pairs=snakemake.params[0],
    )
