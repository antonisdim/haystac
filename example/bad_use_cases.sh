#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

#### haystac config ####

# bad base mismatch base probability
haystac config --mismatch-probability 0.5

# bad scaling value
haystac config --bowtie2-scaling 150

# bad boolean
haystac config --use-conda tr


#### haystac database ####

# create an empty file for some examples
touch empty.txt

# bad database mode
haystac database --mode database --output db_example

# no query whatsoever
haystac database --mode build --output db_example

# bad query
haystac database --mode build --query qwerty --output db_example

# empty query file
haystac database --mode build --query-file empty.txt --output db_example

# empty acessions file
haystac database --mode build --accessions-file empty.txt --output db_example

# bad delimiter accessions file
echo "Yersinia_aldovae  NZ_CP009781.1" > acc_bad_delimiter.txt
haystac database --mode build --accessions-file acc_bad_delimiter.txt --output db_example

# duplicated accession
echo -e "Yersinia_aldovae\tNZ_CP009781.1" > acc_duplicated.txt
echo -e "Yersinia_aldovae\tNZ_CP009781.1" >> acc_duplicated.txt
haystac database --mode build --accessions-file acc_duplicated.txt --output db_example

# query and query file together
echo "Yersinia" > q_file.txt
haystac database --mode build --query "Yersinia" --query-file q_file.txt --output db_example

# refseq rep and mtDNA flags together
haystac database --mode build --refseq-rep --mtDNA --output db_example

# empty exclude-accessions
haystac database --mode build --refseq-rep --output db_example --exclude-accessions

# trying to set an option on top of the other (in this case --refseq-rep)
haystac database --mode build --refseq-rep --output db_example


#### haystac sample ####

# produce some inputs
echo "@ERR1018966.1 HWI-D00554:43:C5LFJANXX:7:1101:1879:2212 length=37
GATCGGCCGCCTCGCGCTCGTCGATCGCCAGGTAGCG
+ERR1018966.1 HWI-D00554:43:C5LFJANXX:7:1101:1879:2212 length=37
B@BCBGG>GGGGGGGGGGGGGGGGGGGGGGBGGEGGC
" > example_input.fastq

cp example_input.fastq example_input_r1.fastq
cp example_input.fastq example_input_r2.fastq

echo ">example_seq1
GATCGGCCGCCTCGCGCTCGTCGATCGCCAGGTAGCG" > example.fasta

# emtpy fastq fill
rm -rf example_sample
haystac sample --output example_sample --fastq empty.txt

# bad path for fastq file
rm -rf example_sample
haystac sample --output example_sample --fastq bad_path.fastq

# provide fasta instead of fastq
rm -rf example_sample
haystac sample --output example_sample --fastq example.fasta

# not providing any seqeunce data for a sample
rm -rf example_sample
haystac sample --output example_sample

# providing only one mate when PE
rm -rf example_sample
haystac sample --output example_sample --fastq-r1 example_input_r1.fastq
haystac sample --output example_sample --fastq-r2 example_input_r2.fastq

# providing inputs for both --input and --input_r1
rm -rf example_sample
haystac sample --output example_sample --fastq example_input.fastq \
  --fastq-r2 example_input_r2.fastq

# providing inputs for both --sra and --fastq
rm -rf example_sample
haystac sample --output example_sample --fastq example_input.fastq --sra ERR1018966

# using --collapse with single end
rm -rf example_sample
haystac sample --output example_sample --fastq example_input.fastq --collapse True

# using --collapse and --trim-adapters False
rm -rf example_sample
haystac sample --output example_sample --fastq-r1 example_input_r1.fastq --fastq-r2 example_input_r2.fastq \
  --collapse True --trim-adapters False

#### haystac analyse ####

# bad mode
haystac analyse --mode likelihood

# not all the required arguments (--mode, --database, --sample, --output)
haystac analyse --mode abundances --output analysis_example

# no db yaml file in db folder
haystac analyse --mode abundances --output analysis_example --database db_example_empty --sample example_sample

# no sample yaml file in sample folder
haystac analyse --mode abundances --output analysis_example --database db_example --sample example_sample_empty