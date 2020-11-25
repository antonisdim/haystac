#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

#### haystac config ####

# bad base mismatch base probability
haystack config --mismatch-probability 0.5

# bad scaling value
haystack config --bowtie2-scaling 150

# bad boolean
haystack config --use-conda tr


#### haystac database ####

# create an empty file for some examples
touch empty.txt

# bad database mode
haystack database --mode database --output db_example

# no query whatsoever
haystack database --mode build --output db_example

# bad query
haystack database --mode build --query qwerty --output db_example

# empty query file
haystack database --mode build --query-file empty.txt --output db_example

# empty acessions file
haystack database --mode build --accessions-file empty.txt --output db_example

# bad delimiter accessions file
echo "Yersinia_aldovae  NZ_CP009781.1" > acc_bad_delimiter.txt
haystack database --mode build --accessions-file acc_bad_delimiter.txt --output db_example

# duplicated accession
echo -e "Yersinia_aldovae\tNZ_CP009781.1" > acc_duplicated.txt
echo -e "Yersinia_aldovae\tNZ_CP009781.1" >> acc_duplicated.txt
haystack database --mode build --accessions-file acc_duplicated.txt --output db_example

# query and query file together
echo "Yersinia" > q_file.txt
haystack database --mode build --query "Yersinia" --query-file q_file.txt --output db_example

# refseq rep and mtDNA flags together
haystack database --mode build --refseq-rep --mtDNA --output db_example

# empty exclude-accessions
haystack database --mode build --refseq-rep --output db_example --exclude-accessions

# trying to set an option on top of the other (in this case --refseq-rep)
haystack database --mode build --refseq-rep --output db_example


#### haystac sample ####

# emtpy fastq fill
rm -rf example_sample
haystack sample --sample-prefix bad_example --output example_sample --fastq empty.txt

# bad path for fastq file
rm -rf example_sample
haystack sample --sample-prefix bad_example --output example_sample --fastq bad_path.fastq

# provide fasta instead of fastq
rm -rf example_sample
haystack sample --sample-prefix bad_example --output example_sample --fastq example.fasta

# not providing any seqeunce data for a sample
rm -rf example_sample
haystack sample --sample-prefix bad_example --output example_sample

# providing only one mate when PE
rm -rf example_sample
haystack sample --sample-prefix bad_example --output example_sample --fastq-r1 example_input_r1.fastq
haystack sample --sample-prefix bad_example --output example_sample --fastq-r2 example_input_r2.fastq

# providing inputs for both --input and --input_r1
rm -rf example_sample
haystack sample --sample-prefix bad_example --output example_sample --fastq example_input.fastq \
  --fastq-r2 example_input_r2.fastq

# providing inputs for both --sra and --fastq
rm -rf example_sample
haystack sample --output example_sample --fastq example_input.fastq --sra ERR1018966

# providing both an --sra and --sample-prefix value
rm -rf example_sample
haystack sample --output example_sample --sra ERR1018966 --sample-prefix bad_example

# using --collapse with single end
rm -rf example_sample
haystack sample --output example_sample --sample-prefix bad_example --fastq example_input.fastq --collapse True

# using --collapse and --trim-adapters False
rm -rf example_sample
haystack sample --output example_sample --fastq-r1 example_input_r1.fastq --fastq-r2 example_input_r2.fastq \
  --sample-prefix bad_example --collapse True --trim-adapters False

# not providing a --sample-prefix
rm -rf example_sample
haystack sample --output example_sample --fastq example_input.fastq





