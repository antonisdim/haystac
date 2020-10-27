#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import random
import sys


def random_db_paths(input_list, output, seed):
    """Function that randomises the genomes that are going to be split into chunks for the bt2 index."""

    random.Random(seed).shuffle(input_list)

    with open(output, "w") as fout:
        for path in input_list:
            fout.write(f"{path}\n")


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    random_db_paths(
        input_list=snakemake.input,
        output=snakemake.output[0],
        seed=snakemake.params[0],
    )
