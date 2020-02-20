#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json

def do_something(input_file, output_file):
    """
    Do something
    """
    
    blah = {}

    with open(output_file, 'w') as fout:
        json.dump(blah, fout)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')
    
    do_something(
		input_file=snakemake.input[0],
        output_file=snakemake.output[0]
    )
