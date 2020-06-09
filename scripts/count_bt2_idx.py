#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import gzip
import os
import sys

from Bio import SeqIO
from Bio import bgzf
from psutil import virtual_memory

MAX_MEM_MB = virtual_memory().total / (1024 ** 2)


# if not MEM_RESOURCES_MB then MAX_MEM_MB

def count_bt2_idx(mem_resources, mem_rescaling_factor, input_file_list, query, output_dir, out_chunk_paths_file):

    mem_resources_mb = mem_resources

    if not mem_resources_mb:
        mem_resources_mb = MAX_MEM_MB

    # read all the files necessary - 1 or 2 - and read the size of file the input function  is returning

    file_size = 0

    input_file_list = [i for i in str(input_file_list).split()]
    for input_file in input_file_list:
        file_size += round(os.stat(input_file).st_size / (1024 ** 2))

    # see if mem_resources are x2.5 its size or based on any other rescaling factor
    scaling_factor = mem_rescaling_factor

    chunk_num = 1
    output_file_path = output_dir + query + "_chunk" + str(chunk_num) + ".fasta.gz"

    output_file_list = [output_file_path]

    if not mem_resources_mb / float(file_size) >= float(scaling_factor):

        fout = bgzf.open(output_file_path, 'wt')
        for input_file in input_file_list:
            binary_handle = bgzf.open(input_file, 'rt')

            for seq_record in SeqIO.parse(binary_handle, "fasta"):
                fout.write(">" + str(seq_record.description) + "\n" + str(seq_record.seq) + "\n")
                size = round(os.stat(output_file_path).st_size / (1024 ** 2))

                if float(size) >= file_size / float(scaling_factor):
                    print("Chunk number {num} has reached its size limit".format(num=chunk_num), file=sys.stderr)
                    fout.close()
                    chunk_num += 1
                    output_file_path = output_dir + query + "_chunk" + str(chunk_num) + ".fasta.gz"
                    output_file_list.append(output_file_path)
                    fout = open(output_file_path, "wt")

    else:

        fout = bgzf.open(output_file_path, 'wt')
        for input_file in input_file_list:
            binary_handle = gzip.open(input_file, 'rt')

            for seq_record in SeqIO.parse(binary_handle, "fasta"):
                fout.write(">" + str(seq_record.description) + "\n" + str(seq_record.seq) + "\n")

    # write the output paths into a file that can be read later

    with open(out_chunk_paths_file, 'w') as outfile:
        print('\n'.join(output_file_list), file=outfile)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    count_bt2_idx(
        mem_resources=snakemake.params[2],
        mem_rescaling_factor=snakemake.params[3],
        input_file_list=snakemake.input,
        query=snakemake.params[0],
        output_dir=snakemake.params[1],
        out_chunk_paths_file=snakemake.output[0],
    )
