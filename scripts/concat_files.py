#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import sys


def concat_files(list_of_files, output_file):
    print("We're concatenating the files ...", file=sys.stderr)

    with open(output_file, 'w') as fout:
        for file in list_of_files:
            print(file, file=sys.stderr)
            with open(file, 'r') as fin:
                shutil.copyfileobj(fin, fout)

    print("{output} file created".format(output=output_file), file=sys.stderr)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    concat_files(
        list_of_files=snakemake.input,
        output_file=snakemake.output[0]
    )
