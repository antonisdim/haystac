#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
from math import ceil

from psutil import virtual_memory

SUBSAMPLE_FIXED_READS = 200000


MEGABYTE = float(1024 ** 2)
MAX_MEM_MB = virtual_memory().total / MEGABYTE


##### Target rules #####

from scripts.rip_utilities import get_total_paths, normalise_name


def get_total_fasta_paths(wildcards):
    """
    Get all the individual fasta file paths for the taxa in our database.
    """

    sequences = get_total_paths(
        wildcards,
        checkpoints,
        config["with_entrez_query"],
        config["refseq_rep"],
        config["with_custom_sequences"],
        config["with_custom_accessions"],
        config["genera"],
    )

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["GBSeq_accession-version"],
        )

        inputs.append(
            config["genome_cache_folder"]
            + "/{orgname}/{accession}.fasta.gz".format(
                orgname=orgname, accession=accession,
            )
        )

    return inputs


checkpoint calculate_bt2_idx_chunks:
    input:
        get_total_fasta_paths,
    log:
        config["db_output"] + "/bowtie/bt2_idx_chunk_num.log",
    output:
        config["db_output"] + "/bowtie/bt2_idx_chunk_num.txt",
    params:
        query=config["db_output"],
        mem_resources=float(config["mem"]),
        mem_rescale_factor=config["bowtie2_scaling"],
    message:
        "The number of index chunks for the filtering alignment are being calculated for query {params.query}. "
        "The size rescaling factor for the chunk is {params.mem_rescale_factor} for the given memory "
        "resources {params.mem_resources}. The log file can be found in {log}."
    script:
        "../scripts/calculate_bt2_idx_chunks.py"


def get_bt2_idx_filter_chunk(wildcards):

    """Pick the files for the specific bt2 index chunk"""

    chunk_files = []
    chunk_size = float(config["mem"]) / float(config["bowtie2_scaling"])
    total_size = 0.0

    for fasta_file in get_total_fasta_paths(wildcards):
        total_size += os.stat(fasta_file).st_size / MEGABYTE

        chunk = (total_size // chunk_size) + 1

        if chunk == int(wildcards.chunk_num):
            chunk_files.append(fasta_file)
        elif chunk > int(wildcards.chunk_num):
            break

    return chunk_files


rule create_bt2_idx_filter_chunk:
    input:
        get_bt2_idx_filter_chunk,
    log:
        config["db_output"] + "/bowtie/bt2_idx_filter_{chunk_num}.log",
    output:
        config["db_output"] + "/bowtie/chunk{chunk_num}.fasta.gz",
    message:
        "Creating fasta chunk {wildcards.chunk_num} for the filtering alignment bowtie2 index for "
        "the given query.The output can be found in {output} and its log can be found in {log}."
    script:
        "../scripts/bowtie2_multifasta.py"


rule bowtie_index:
    input:
        fasta_chunk=config["db_output"] + "/bowtie/chunk{chunk_num}.fasta.gz",
    log:
        config["db_output"] + "/bowtie/{chunk_num}_index.log",
    output:
        expand(
            config["db_output"] + "/bowtie/chunk{chunk_num}.{n}.bt2l",
            n=[1, 2, 3, 4],
            allow_missing=True,
        ),
        expand(
            config["db_output"] + "/bowtie/chunk{chunk_num}.rev.{n}.bt2l",
            n=[1, 2],
            allow_missing=True,
        ),
    benchmark:
        repeat("benchmarks/bowtie_index_chunk{chunk_num}.benchmark.txt", 1)
    message:
        "Bowtie2 index for chunk {input.fasta_chunk} is being built. The log file can be found in {log}."
    shell:
        "bowtie2-build --large-index {input.fasta_chunk} "+
         config['db_output']+"/bowtie/chunk{wildcards.chunk_num} &> {log}"


def get__bt2_idx_chunk_paths(wildcards):

    """Get the paths for the index chunks for the filtering bowtie2 alignment"""

    get_chunk_num = checkpoints.calculate_bt2_idx_chunks.get()
    idx_chunk_total = ceil(float(open(get_chunk_num.output[0]).read().strip()))

    return expand(
        config["db_output"] + "/bowtie/chunk{chunk_num}.1.bt2l",
        chunk_num=[x + 1 if idx_chunk_total > 1 else 1 for x in range(idx_chunk_total)],
    )


rule bowtie_index_done:
    input:
        get__bt2_idx_chunk_paths,
    output:
        config["db_output"] + "/bowtie/bowtie_index.done",
    benchmark:
        repeat("benchmarks/bowtie_index_done", 1)
    message:
        "The bowtie2 indices for all the chunks {input} have been built."
    shell:
        "touch {output}"


def get_inputs_for_bowtie_r1(wildcards):

    if config["trim_adapters"]:
        if config["PE_MODERN"]:
            return config[
                "sample_output_dir"
            ] + "/fastq_inputs/PE_mod/{sample}_R1_adRm.fastq.gz".format(
                sample=wildcards.sample
            )
        elif config["PE_ANCIENT"]:
            return config[
                "sample_output_dir"
            ] + "/fastq_inputs/PE_anc/{sample}_adRm.fastq.gz".format(
                sample=wildcards.sample
            )
        elif config["SE"]:
            return config[
                "sample_output_dir"
            ] + "/fastq_inputs/SE/{sample}_adRm.fastq.gz".format(
                sample=wildcards.sample
            )

    if config["PE_MODERN"]:
        return config["fastq_R1"]
    elif config["PE_ANCIENT"] or config["SE"]:
        return config["fastq"]


def get_inputs_for_bowtie_r2(wildcards):

    if config["trim_adapters"]:
        if config["PE_MODERN"]:
            return "fastq_inputs/PE_mod/{sample}_R2_adRm.fastq.gz".format(
                sample=wildcards.sample
            )

    if config["PE_MODERN"]:
        return config["fastq_R2"]


rule bowtie_alignment_single_end:
    input:
        fastq=get_inputs_for_bowtie_r1,
        bt2idx=config["db_output"] + "/bowtie/chunk{chunk_num}.1.bt2l",
    log:
        config["analysis_output_dir"] + "/bam/{sample}_chunk{chunk_num}.log",
    output:
        bam_file=temp(
            config["analysis_output_dir"] + "/bam/SE_{sample}_sorted_chunk{chunk_num}.bam"
        ),
    benchmark:
        repeat("benchmarks/bowtie_alignment_{sample}_chunk{chunk_num}.benchmark.txt", 1)
    params:
        index=config["db_output"] + "/bowtie/chunk{chunk_num}",
    threads: config["bowtie2_threads"]
    message:
        "The filtering alignment for file {input.fastq}, of sample {wildcards.sample}, "
        "for index chunk number {wildcards.chunk_num} is being executed, "
        "with number {threads} of threads. The log file can be found in {log}."
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -U {input.fastq} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"


rule bowtie_alignment_paired_end:
    input:
        fastq_r1=get_inputs_for_bowtie_r1,
        fastq_r2=get_inputs_for_bowtie_r2,
        bt2idx=config["db_output"] + "/bowtie/chunk{chunk_num}.1.bt2l",
    log:
        config["analysis_output_dir"] + "/bam/{sample}_chunk{chunk_num}.log",
    output:
        bam_file=temp(
            config["analysis_output_dir"] + "/bam/PE_{sample}_sorted_chunk{chunk_num}.bam"
        ),
    benchmark:
        repeat("benchmarks/bowtie_alignment_{sample}_chunk{chunk_num}.benchmark.txt", 1)
    params:
        index=config["db_output"] + "/bowtie/chunk{chunk_num}",
    threads: config["bowtie2_threads"]
    message:
        "The filtering alignment for files {input.fastq_r1} and {input.fastq_r2}, of sample {wildcards.sample}, "
        "for index chunk number {wildcards.chunk_num} is being executed, "
        "with number {threads} of threads. The log file can be found in {log}."
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -1 {input.fastq_r1} -2 {input.fastq_r2} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"


def get_sorted_bam_paths(wildcards):
    """Get the sorted bam paths for each index chunk"""

    # get_chunk_num = checkpoints.calculate_bt2_idx_chunks.get(query=wildcards.query)
    # idx_chunk_total = ceil(float(open(get_chunk_num.output[0]).read().strip()))
    idx_chunk_total = ceil(
        float(
            open(config["db_output"] + "/bowtie/bt2_idx_chunk_num.txt").read().strip()
        )
    )

    reads = ["PE"] if config["PE_MODERN"] else ["SE"]

    return expand(
        config["analysis_output_dir"]
        + "/bam/{reads}_{sample}_sorted_chunk{chunk_num}.bam",
        reads=reads,
        sample=wildcards.sample,
        chunk_num=range(1, idx_chunk_total + 1),
    )


rule merge_bams:
    input:
        aln_path=get_sorted_bam_paths,
    log:
        config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_merge_bams.log",
    output:
        config["analysis_output_dir"] + "/bam/{read_mode}_{sample}_sorted.bam",
    message:
        "Merging the bam files ({input}) produced by the filtering alignment stage for sample {wildcards.sample}. "
        "The log file can be found in {log}."
    shell:
        "samtools merge -f {output} {input.aln_path} 2> {log}"


rule extract_fastq_single_end:
    input:
        config["analysis_output_dir"] + "/bam/SE_{sample}_sorted.bam",
    log:
        config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.log",
    output:
        config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.fastq.gz",
    benchmark:
        repeat("benchmarks/extract_fastq_single_end_{sample}.benchmark.txt", 1)
    message:
        "Extracting all the aligned reads for sample {wildcards.sample} and storing them in {output}. "
        "The log file can be found in {log}."
    shell:
        "( samtools view -h -F 4 {input} | samtools fastq -c 6 - | seqkit rmdup -n -o {output} ) 2> {log}"


rule extract_fastq_paired_end:
    input:
        config["analysis_output_dir"] + "/bam/PE_{sample}_sorted.bam",
    log:
        config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq.log",
    output:
        config["analysis_output_dir"] + "/fastq/PE/{sample}_R1_mapq.fastq.gz",
        config["analysis_output_dir"] + "/fastq/PE/{sample}_R2_mapq.fastq.gz",
    benchmark:
        repeat("benchmarks/extract_fastq_paired_end_{sample}.benchmark.txt", 1)
    message:
        "Extracting all the aligned reads for sample {wildcards.sample} and storing them in {output}. "
        "The log file can be found in {log}."
    shell:
        "( samtools view -h -F 4 {input} "
        "| samtools fastq -c 6 -1 " + config[
            "analysis_output_dir"
        ] + "/fastq/PE/{wildcards.sample}_temp_R1.fastq.gz "
        "-2 " + config[
            "analysis_output_dir"
        ] + "/fastq/PE/{wildcards.sample}_temp_R2.fastq.gz -0 /dev/null -s /dev/null -;"
        "seqkit rmdup -n " + config[
            "analysis_output_dir"
        ] + "/fastq/PE/{wildcards.sample}_temp_R1.fastq.gz -o {output[0]}; "
        "seqkit rmdup -n " + config[
            "analysis_output_dir"
        ] + "/fastq/PE/{wildcards.sample}_temp_R2.fastq.gz -o {output[1]}; "
        "rm " + config[
            "analysis_output_dir"
        ] + "/fastq/PE/{wildcards.sample}_temp_R1.fastq.gz; "
        "rm " + config[
            "analysis_output_dir"
        ] + "/fastq/PE/{wildcards.sample}_temp_R2.fastq.gz ) 2> {log}"


rule average_fastq_read_len_single_end:
    input:
        config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.fastq.gz",
    log:
        config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq_readlen.log",
    output:
        config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.readlen",
    benchmark:
        repeat("benchmarks/average_fastq_read_len_single_end_{sample}.benchmark.txt", 1)
    message:
        "Calculating the average read length for sample {wildcards.sample} from file {input} "
        "and storing its value in {output}. The log file can be found in {log}."
    shell:
        "( seqtk sample {input} {SUBSAMPLE_FIXED_READS} | seqtk seq -A | grep -v '^>'"
        "| awk '{{count++; bases += length}} END {{print bases/count}}' 1> {output} ) 2> {log}"


rule average_fastq_read_len_paired_end:
    input:
        mate1=config["analysis_output_dir"] + "/fastq/PE/{sample}_R1_mapq.fastq.gz",
        mate2=config["analysis_output_dir"] + "/fastq/PE/{sample}_R2_mapq.fastq.gz",
    log:
        config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq_readlen.log",
    output:
        mate1=temp(config["analysis_output_dir"] + "/fastq/{sample}_R1_mapq.readlen"),
        mate2=temp(config["analysis_output_dir"] + "/fastq/{sample}_R2_mapq.readlen"),
        pair=config["analysis_output_dir"] + "/fastq/PE/{sample}_mapq_pair.readlen",
    benchmark:
        repeat("benchmarks/average_fastq_read_len_paired_end_{sample}.benchmark.txt", 1)
    message:
        "Calculating the average read length for sample {wildcards.sample} from files {input.mate1} and {input.mate2} "
        "and storing its value in {output}. The log file can be found in {log}."
    shell:
        "( seqtk sample {input.mate1} {SUBSAMPLE_FIXED_READS} | seqtk seq -A | grep -v '^>' "
        "| awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output.mate1} ) 2>> {log}; "
        "( seqtk sample {input.mate2} {SUBSAMPLE_FIXED_READS} | seqtk seq -A | grep -v '^>' "
        "| awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output.mate2} ) 2>> {log}; "
        "( cat {output.mate1} {output.mate2} "
        "| awk '{{sum += $1; n++ }} END {{if (n>0) print sum/n;}}' 1> {output.pair} ) 2>> {log} "
