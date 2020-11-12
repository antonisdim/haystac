#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import shutil
import sys


# TODO get rid of this file and do the concat operations on the command line using `cat {list} | xargs cat > {output}`
def concat_files(list_of_files, output_file):
    """Concatenating files that are provided as input into a single file."""

    print("Concatenating files ...", file=sys.stderr)

    with open(output_file, "w") as fout:
        for file in list_of_files:
            print(file, file=sys.stderr)
            with open(file, "r") as fin:
                shutil.copyfileobj(fin, fout)

    print(f"{output_file} file created", file=sys.stderr)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    concat_files(list_of_files=snakemake.input, output_file=snakemake.output[0])
