#!/usr/bin/env python
# -*- coding: utf-8 -*-

from multiprocessing import cpu_count
from psutil import virtual_memory
import pandas as pd

SUBSAMPLE_FIXED_READS = 200000
WITH_REFSEQ_REP = config["WITH_REFSEQ_REP"]
WITH_ENTREZ_QUERY = config["WITH_ENTREZ_QUERY"]
WITH_CUSTOM_SEQUENCES = config["WITH_CUSTOM_SEQUENCES"]
WITH_CUSTOM_ACCESSIONS = config["WITH_CUSTOM_ACCESSIONS"]
SRA_LOOKUP = config["SRA_LOOKUP"]
PE_ANCIENT = config["PE_ANCIENT"]
PE_MODERN = config["PE_MODERN"]
SE = config["SE"]
MAX_MEM_MB = virtual_memory().total / (1024 ** 2)
MEM_RESOURCES_MB = float(config["MEM_RESOURCES_MB"])
MEM_RESCALING_FACTOR = config["MEM_RESCALING_FACTOR"]
WITH_DATA_PREPROCESSING = config["WITH_DATA_PREPROCESSING"]

##### Target rules #####


def get_total_fasta_paths(wildcards):
    """
    Get all the individual fasta file paths for the taxa in our database.
    """

    sequences = pd.DataFrame()

    if WITH_ENTREZ_QUERY:
        pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
        sequences = pd.read_csv(pick_sequences.output[0], sep="\t")

        if len(sequences) == 0:
            raise RuntimeError("The entrez pick sequences file is empty.")

    if WITH_REFSEQ_REP:
        refseq_rep_prok = checkpoints.entrez_refseq_accessions.get(
            query=wildcards.query
        )

        refseq_genomes = pd.read_csv(refseq_rep_prok.output[0], sep="\t")
        genbank_genomes = pd.read_csv(refseq_rep_prok.output[1], sep="\t")
        assemblies = pd.read_csv(refseq_rep_prok.output[2], sep="\t")
        refseq_plasmids = pd.read_csv(refseq_rep_prok.output[3], sep="\t")
        genbank_plasmids = pd.read_csv(refseq_rep_prok.output[4], sep="\t")

        invalid_assemblies = checkpoints.entrez_invalid_assemblies.get(
            query=wildcards.query
        )
        invalid_assembly_sequences = pd.read_csv(invalid_assemblies.output[0], sep="\t")

        assemblies = assemblies[
            ~assemblies["GBSeq_accession-version"].isin(
                invalid_assembly_sequences["GBSeq_accession-version"]
            )
        ]

        if WITH_ENTREZ_QUERY:
            sequences = pd.concat(
                [
                    sequences,
                    refseq_genomes,
                    genbank_genomes,
                    assemblies,
                    refseq_plasmids,
                    genbank_plasmids,
                ]
            )
        else:
            sequences = pd.concat(
                [
                    refseq_genomes,
                    genbank_genomes,
                    assemblies,
                    refseq_plasmids,
                    genbank_plasmids,
                ]
            )

    if WITH_CUSTOM_SEQUENCES:
        custom_fasta_paths = pd.read_csv(
            config["custom_seq_file"],
            sep="\t",
            header=None,
            names=["species", "GBSeq_accession-version", "path"],
        )

        custom_seqs = custom_fasta_paths[["species", "GBSeq_accession-version"]]

        sequences = sequences.append(custom_seqs)

    if WITH_CUSTOM_ACCESSIONS:
        custom_accessions = pd.read_csv(
            config["custom_acc_file"],
            sep="\t",
            header=None,
            names=["species", "GBSeq_accession-version"],
        )

        sequences = sequences.append(custom_accessions)

    inputs = []

    if SPECIFIC_GENUS:
        sequences = sequences[
            sequences["species"].str.contains("|".join(SPECIFIC_GENUS))
        ]

    for key, seq in sequences.iterrows():
        orgname, accession = (
            seq["species"].replace(" ", "_").replace("[", "").replace("]", ""),
            seq["GBSeq_accession-version"],
        )

        inputs.append(
            "database/{orgname}/{accession}.fasta.gz".format(
                orgname=orgname, accession=accession,
            )
        )

    return inputs


checkpoint calculate_bt2_idx_chunks:
    input:
        get_total_fasta_paths,
    log:
        "{query}/bowtie/{query}_bt2_idx_chunk_num.log",
    output:
        "{query}/bowtie/{query}_bt2_idx_chunk_num.txt",
    params:
        query="{query}",
        mem_resources=MEM_RESOURCES_MB,
        mem_rescaling_factor=MEM_RESCALING_FACTOR,
    message:
        "The number of index chunks for the filtering alignment are being calculated for query {params.query}. "
        "The size rescaling factor for the chunk is {params.mem_rescaling_factor} for the given memory "
        "resources {params.mem_resources}. The log file can be found in {log}."
    script:
        "../scripts/calculate_bt2_idx_chunks.py"


def get_bt2_idx_filter_chunk(wildcards):

    """Pick the files for the specific bt2 index chunk"""

    chunk = 1
    chunk_to_chose = int(wildcards.chunk_num)

    mem_resources_mb = MEM_RESOURCES_MB

    if not mem_resources_mb:
        mem_resources_mb = MAX_MEM_MB

    incremental_file_size = float(0)

    input_file_list = get_total_fasta_paths(wildcards)

    chunk_path_list = []

    for input_file in input_file_list:

        size_to_be = incremental_file_size + os.stat(input_file).st_size / float(
            1024 ** 2
        )

        if incremental_file_size == 0:
            chunk_path_list = []
            incremental_file_size += os.stat(input_file).st_size / float(1024 ** 2)
            chunk_path_list.append(input_file)

        elif size_to_be <= mem_resources_mb / float(MEM_RESCALING_FACTOR):
            incremental_file_size += os.stat(input_file).st_size / float(1024 ** 2)
            chunk_path_list.append(input_file)

        elif size_to_be > mem_resources_mb / float(MEM_RESCALING_FACTOR):
            if chunk == chunk_to_chose:
                break
            incremental_file_size = 0
            chunk += 1
            incremental_file_size += os.stat(input_file).st_size / float(1024 ** 2)
            chunk_path_list = [input_file]

    return chunk_path_list


rule create_bt2_idx_filter_chunk:
    input:
        get_bt2_idx_filter_chunk,
    log:
        "{query}/bowtie/{query}_bt2_idx_filter_{chunk_num}.log",
    output:
        "{query}/bowtie/{query}_chunk{chunk_num}.fasta.gz",
    message:
        "Creating fasta chunk {wildcards.chunk_num} for the filtering alignment bowtie2 index for "
        "query {wildcards.query}.The output can be found in {output} and its log can be found in {log}."
    script:
        "../scripts/bowtie2_multifasta.py"


rule bowtie_index:
    input:
        fasta_chunk="{query}/bowtie/{query}_chunk{chunk_num}.fasta.gz",
    log:
        "{query}/bowtie/{query}_{chunk_num}_index.log",
    output:
        expand("{{query}}/bowtie/{{query}}_chunk{{chunk_num}}.{n}.bt2l", n=[1, 2, 3, 4]),
        expand("{{query}}/bowtie/{{query}}_chunk{{chunk_num}}.rev.{n}.bt2l", n=[1, 2]),
    benchmark:
        repeat("benchmarks/bowtie_index_{query}_chunk{chunk_num}.benchmark.txt", 1)
    message:
        "Bowtie2 index for chunk {input.fasta_chunk} is being built. The log file can be found in {log}."
    shell:
        "bowtie2-build --large-index {input.fasta_chunk} "
        "{wildcards.query}/bowtie/{wildcards.query}_chunk{wildcards.chunk_num} &> {log}"


def get__bt2_idx_chunk_paths(wildcards):

    get_chunk_num = checkpoints.calculate_bt2_idx_chunks.get(query=wildcards.query)
    idx_chunk_total = int(float(open(get_chunk_num.output[0]).read()))

    return expand(
        "{query}/bowtie/{query}_chunk{chunk_num}.1.bt2l",
        query=wildcards.query,
        chunk_num=[x + 1 if idx_chunk_total > 1 else 1 for x in range(idx_chunk_total)],
    )


rule bowtie_index_done:
    input:
        get__bt2_idx_chunk_paths,
    log:
        "{query}/bowtie/bowtie_index.log",
    output:
        "{query}/bowtie/bowtie_index.done",
    benchmark:
        repeat("benchmarks/bowtie_index_done_{query}", 1)
    message:
        "The bowtie2 indices for all the chunks {input} have been built. The log file can be found in {log}."
    shell:
        "touch {output} 2> {log}"


def get_inputs_for_bowtie_r1(wildcards):
    print(wildcards.sample)

    if SRA_LOOKUP:
        if PE_MODERN:
            return "fastq_inputs/PE_mod/{sample}_R1_adRm.fastq.gz".format(
                sample=wildcards.sample
            )
        elif PE_ANCIENT:
            return "fastq_inputs/PE_anc/{sample}_adRm.fastq.gz".format(
                sample=wildcards.sample
            )
        elif SE:
            return "fastq_inputs/SE/{sample}_adRm.fastq.gz".format(
                sample=wildcards.sample
            )

    else:
        if PE_MODERN:
            return config["sample_fastq_R1"]
        elif PE_ANCIENT:
            if WITH_DATA_PREPROCESSING:
                return "fastq_inputs/PE_anc/{sample}_adRm.fastq.gz".format(
                    sample=wildcards.sample
                )
            else:
                return config["sample_fastq"]
        elif SE:
            if WITH_DATA_PREPROCESSING:
                return "fastq_inputs/SE/{sample}_adRm.fastq.gz".format(
                    sample=wildcards.sample
                )
            else:
                return config["sample_fastq"]


def get_inputs_for_bowtie_r2(wildcards):
    if SRA_LOOKUP:
        if PE_MODERN:
            return "fastq_inputs/PE_mod/{sample}_R2_adRm.fastq.gz".format(
                sample=wildcards.sample
            )
        else:
            return ""

    else:
        if PE_MODERN:
            return config["sample_fastq_R2"]
        else:
            return ""


rule bowtie_alignment_single_end:
    input:
        fastq=get_inputs_for_bowtie_r1,
        bt2idx="{query}/bowtie/{query}_chunk{chunk_num}.1.bt2l",
    log:
        "{query}/bam/{sample}_chunk{chunk_num}.log",
    output:
        bam_file=temp("{query}/bam/SE_{sample}_sorted_chunk{chunk_num}.bam"),
    benchmark:
        repeat("benchmarks/bowtie_alignment_{query}_{sample}_chunk{chunk_num}.benchmark.txt", 1)
    params:
        index="{query}/bowtie/{query}_chunk{chunk_num}",
    threads: config["bowtie2_treads"]
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
        bt2idx="{query}/bowtie/{query}_chunk{chunk_num}.1.bt2l",
    log:
        "{query}/bam/{sample}_chunk{chunk_num}.log",
    output:
        bam_file=temp("{query}/bam/PE_{sample}_sorted_chunk{chunk_num}.bam"),
    benchmark:
        repeat("benchmarks/bowtie_alignment_{query}_{sample}_chunk{chunk_num}.benchmark.txt", 1)
    params:
        index="{query}/bowtie/{query}_chunk{chunk_num}",
    threads: config["bowtie2_treads"]
    message:
        "The filtering alignment for files {input.fastq_r1} and {input.fastq_r2}, of sample {wildcards.sample}, "
        "for index chunk number {wildcards.chunk_num} is being executed, "
        "with number {threads} of threads. The log file can be found in {log}."
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -1 {input.fastq_r1} -2 {input.fastq_r2} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"


def get_sorted_bam_paths(wildcards):

    get_chunk_num = checkpoints.calculate_bt2_idx_chunks.get(query=wildcards.query)
    idx_chunk_total = int(open(get_chunk_num.output[0]).read())

    if PE_MODERN:
        return expand(
            "{query}/bam/{reads}_{sample}_sorted_chunk{chunk_num}.bam",
            query=wildcards.query,
            reads=["PE"],
            sample=wildcards.sample,
            chunk_num=[
                x + 1 if idx_chunk_total > 1 else 1 for x in range(idx_chunk_total)
            ],
        )
    if PE_ANCIENT or SE:
        return expand(
            "{query}/bam/{reads}_{sample}_sorted_chunk{chunk_num}.bam",
            query=wildcards.query,
            reads=["SE"],
            sample=wildcards.sample,
            chunk_num=[
                x + 1 if idx_chunk_total > 1 else 1 for x in range(idx_chunk_total)
            ],
        )


rule merge_bams_single_end:
    input:
        aln_path=get_sorted_bam_paths, # aln_done="{query}/bam/{sample}_all_aln.done",
    log:
        "{query}/bam/{sample}_merge_bams.log",
    output:
        "{query}/bam/SE_{sample}_sorted.bam",
    message:
        "Merging bam files ({input}) produced by the filtering alignment stage for sample {wildcards.sample}. "
        "The log file can be found in {log}."
    shell:
        "samtools merge -f {output} {input.aln_path} 2> {log}"


rule merge_bams_paired_end:
    input:
        aln_path=get_sorted_bam_paths, # aln_done="{query}/bam/{sample}_all_aln.done",
    log:
        "{query}/bam/{sample}_merge_bams.log",
    output:
        "{query}/bam/PE_{sample}_sorted.bam",
    message:
        "Merging bam files ({input}) produced by the filtering alignment stage for sample {wildcards.sample}. "
        "The log file can be found in {log}."
    shell:
        "samtools merge -f {output} {input.aln_path} 2> {log}"


rule extract_fastq_single_end:
    input:
        "{query}/bam/SE_{sample}_sorted.bam",
    log:
        "{query}/fastq/SE/{sample}_mapq.log",
    output:
        "{query}/fastq/SE/{sample}_mapq.fastq.gz",
    benchmark:
        repeat("benchmarks/extract_fastq_single_end_{query}_{sample}.benchmark.txt", 1)
    message:
        "Extracting all the aligned reads for sample {wildcards.sample} and storing them in {output}. "
        "The log file can be found in {log}."
    shell:
        "( samtools view -h -F 4 {input} | samtools fastq -c 6 - | seqkit rmdup -n -o {output} ) 2> {log}" # "( samtools view -h -F 4 {input} | samtools fastq -c 6 - > {output} ) 2> {log}"


rule extract_fastq_paired_end:
    input:
        "{query}/bam/PE_{sample}_sorted.bam",
    log:
        "{query}/fastq/PE/{sample}_mapq.log",
    output:
        "{query}/fastq/PE/{sample}_R1_mapq.fastq.gz",
        "{query}/fastq/PE/{sample}_R2_mapq.fastq.gz",
    benchmark:
        repeat("benchmarks/extract_fastq_paired_end_{query}_{sample}.benchmark.txt", 1)
    message:
        "Extracting all the aligned reads for sample {wildcards.sample} and storing them in {output}. "
        "The log file can be found in {log}."
    shell:
        "( samtools view -h -F 4 {input} "
        "| samtools fastq -c 6 -1 {wildcards.query}/fastq/PE/temp_R1.fastq.gz "
        "-2 {wildcards.query}/fastq/PE/temp_R2.fastq.gz -0 /dev/null -s /dev/null -;"
        "seqkit rmdup -n {wildcards.query}/fastq/PE/{wildcards.sample}_temp_R1.fastq.gz -o {output[0]}; "
        "seqkit rmdup -n {wildcards.query}/fastq/PE/{wildcards.sample}_temp_R2.fastq.gz -o {output[1]}; "
        "rm {wildcards.query}/fastq/PE/{wildcards.sample}_temp_R1.fastq.gz; "
        "rm {wildcards.query}/fastq/PE/{wildcards.sample}_temp_R2.fastq.gz ) 2> {log}" # "( samtools view -h -F 4 {input} "
         # "| samtools fastq -c 6 -1 {output[0]} -2 {output[1]} -0 /dev/null -s /dev/null - ) 2> {log}"


# ruleorder: extract_fastq_paired_end > extract_fastq_single_end


rule average_fastq_read_len_single_end:
    input:
        "{query}/fastq/SE/{sample}_mapq.fastq.gz",
    log:
        "{query}/fastq/SE/{sample}_mapq_readlen.log",
    output:
        "{query}/fastq/SE/{sample}_mapq.readlen",
    benchmark:
        repeat("benchmarks/average_fastq_read_len_single_end_{query}_{sample}.benchmark.txt", 1)
    params:
        sample_size=SUBSAMPLE_FIXED_READS,
    message:
        "Calculating the average read length for sample {wildcards.sample} from file {input} "
        "and storing its value in {output}. The log file can be found in {log}."
    shell:
        "seqtk sample {input} {params.sample_size} | seqtk seq -A | grep -v '^>' | "
        "awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output} 2> {log}"


rule average_fastq_read_len_paired_end:
    input:
        mate1="{query}/fastq/PE/{sample}_R1_mapq.fastq.gz",
        mate2="{query}/fastq/PE/{sample}_R2_mapq.fastq.gz",
    log:
        "{query}/fastq/PE/{sample}_mapq_readlen.log",
    output:
        mate1=temp("{query}/fastq/{sample}_R1_mapq.readlen"),
        mate2=temp("{query}/fastq/{sample}_R2_mapq.readlen"),
        pair="{query}/fastq/PE/{sample}_mapq_pair.readlen",
    benchmark:
        repeat("benchmarks/average_fastq_read_len_paired_end_{query}_{sample}.benchmark.txt", 1)
    params:
        sample_size=SUBSAMPLE_FIXED_READS,
    message:
        "Calculating the average read length for sample {wildcards.sample} from files {input.mate1} and {input.mate2} "
        "and storing its value in {output}. The log file can be found in {log}."
    shell:
        "seqtk sample {input.mate1} {params.sample_size} | seqtk seq -A | grep -v '^>' | "
        "awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output.mate1} 2> {log}; "
        "seqtk sample {input.mate2} {params.sample_size} | seqtk seq -A | grep -v '^>' | "
        "awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output.mate2} 2> {log}; "
        "cat {output.mate1} {output.mate2} | awk '{{sum += $1; n++ }} END {{if (n>0) print sum/n;}}' "
        "1> {output.pair} 2> {log} "
