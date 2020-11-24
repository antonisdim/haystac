#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

#### haystac config ####

# bash strict mode
set -euo pipefail

# bad base mismatch base probability
haystack config --mismatch-probability 0.5

# bad scaling value
haystack config --bowtie2-scaling 150

# bad boolean
haystack config --use-conda tr


#### haystac database ####

# bad database mode
haystack database --mode database

# no query whatsoever
haystack database --mode build --output db_example

# bad query
haystack database --mode build --query qwerty

# empty query file
touch empty.txt
haystack database --mode build --query-file empty.txt

# empty acessions file
--accessions-file

# bad delimiter accessions file

# duplicated accession



--sequences-file

--exclude-accessions

args.mtDNA and args.refseq_rep

 if args.query and args.query_file:
            raise ValidationError("Please specify either `--query <query>` or `--query-file <path>` but not both.")


