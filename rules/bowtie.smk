#!/usr/bin/env python
# -*- coding: utf-8 -*-

from multiprocessing import cpu_count
from psutil import virtual_memory
import pandas as pd

SUBSAMPLE_FIXED_READS = 200000
WITH_REFSEQ_REP = config['WITH_REFSEQ_REP']
WITH_ENTREZ_QUERY = config['WITH_ENTREZ_QUERY']
SRA_LOOKUP = config['SRA_LOOKUP']
PE_ANCIENT = config['PE_ANCIENT']
PE_MODERN = config['PE_MODERN']
SE = config['SE']
MAX_MEM_MB = virtual_memory().total / (1024 ** 2)
MEM_RESOURCES_MB = config['MEM_RESOURCES_MB']
MEM_RESCALING_FACTOR = config['MEM_RESCALING_FACTOR']
WITH_DATA_PREPROCESSING = config['WITH_DATA_PREPROCESSING']

##### Target rules #####


def filtering_bowtie_aln_inputs(wildcards):
    if WITH_REFSEQ_REP and WITH_ENTREZ_QUERY:
        return ["{query}/bowtie/{query}_entrez.fasta.gz".format(query=wildcards.query),
                "{query}/bowtie/{query}_refseq_prok.fasta.gz".format(query=wildcards.query)]
    elif WITH_ENTREZ_QUERY:
        return ["{query}/bowtie/{query}_entrez.fasta.gz".format(query=wildcards.query)]
    else:
        return ["{query}/bowtie/{query}_refseq_prok.fasta.gz".format(query=wildcards.query)]



checkpoint count_bt2_idx:
    input:
        filtering_bowtie_aln_inputs
    log:
        "{query}/bowtie/{query}_count_bt2_idx.log"
    output:
        "{query}/bowtie/{query}_bt2_idx_chunk_paths.txt"
    params:
        query="{query}",
        output_dir="{query}/bowtie/",
        mem_resources=MEM_RESOURCES_MB,
        mem_rescaling_factor=MEM_RESCALING_FACTOR
    script:
        "../scripts/count_bt2_idx.py"



rule bowtie_index:
    input:
        chunk_path_file="{query}/bowtie/{query}_bt2_idx_chunk_paths.txt",
        fasta_chunk="{query}/bowtie/{query}_chunk{chunk_num}.fasta.gz"
    log:
        "{query}/bowtie/{query}_{chunk_num}_index.log"
    output:
        expand("{{query}}/bowtie/{{query}}_chunk{{chunk_num}}.{n}.bt2l", n=[1, 2, 3, 4]),
        expand("{{query}}/bowtie/{{query}}_chunk{{chunk_num}}.rev.{n}.bt2l", n=[1, 2])
    benchmark:
        repeat("benchmarks/bowtie_index_{query}_chunk{chunk_num}.benchmark.txt", 1)
    # run:
    #     if WITH_REFSEQ_REP and WITH_ENTREZ_QUERY:
    #         shell("cat {input} > {wildcards.query}/bowtie/{wildcards.query}.fasta.gz; "
    #               "bowtie2-build --large-index {wildcards.query}/bowtie/{wildcards.query}.fasta.gz "
    #               "{wildcards.query}/bowtie/{wildcards.query} &> {log}")
    #     else:
    #         shell("bowtie2-build --large-index {input} {wildcards.query}/bowtie/{wildcards.query} &> {log}")
    shell:
        "bowtie2-build --large-index {input.fasta_chunk} "
        "{wildcards.query}/bowtie/{wildcards.query}_chunk{wildcards.chunk_num} &> {log}"



def get_inputs_for_bowtie_r1(wildcards):
    print(wildcards.sample)

    if SRA_LOOKUP:
        if PE_MODERN:
            return "fastq_inputs/PE_mod/{sample}_R1_adRm.fastq.gz".format(sample=wildcards.sample)
        elif PE_ANCIENT:
            return "fastq_inputs/PE_anc/{sample}_adRm.fastq.gz".format(sample=wildcards.sample)
        elif SE:
            return "fastq_inputs/SE/{sample}_adRm.fastq.gz".format(sample=wildcards.sample)

    else:
        if PE_MODERN:
            return config['sample_fastq_R1']
        elif PE_ANCIENT:
            if WITH_DATA_PREPROCESSING:
                return "fastq_inputs/PE_anc/{sample}_adRm.fastq.gz".format(sample=wildcards.sample)
            else:
                return config['sample_fastq']
        elif SE:
            if WITH_DATA_PREPROCESSING:
                return "fastq_inputs/SE/{sample}_adRm.fastq.gz".format(sample=wildcards.sample)
            else:
                return config['sample_fastq']



def get_inputs_for_bowtie_r2(wildcards):
    if SRA_LOOKUP:
        if PE_MODERN:
            return "fastq_inputs/PE_mod/{sample}_R2_adRm.fastq.gz".format(sample=wildcards.sample)
        else:
            return ""

    else:
        if PE_MODERN:
            return config['sample_fastq_R2']
        else:
            return ""



rule bowtie_alignment_single_end:
    input:
        fastq=get_inputs_for_bowtie_r1,
        bt2idx="{query}/bowtie/{query}_chunk{chunk_num}.1.bt2l"
    log:
        "{query}/bam/{sample}_chunk{chunk_num}.log"
    output:
        bam_file=temp("{query}/bam/SE_{sample}_sorted_chunk{chunk_num}.bam")
    benchmark:
        repeat("benchmarks/bowtie_alignment_{query}_{sample}_chunk{chunk_num}.benchmark.txt", 1)
    params:
        index="{query}/bowtie/{query}_chunk{chunk_num}"
    threads:
        cpu_count()
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -U {input.fastq} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"



rule bowtie_alignment_paired_end:
    input:
        fastq_r1=get_inputs_for_bowtie_r1,
        fastq_r2=get_inputs_for_bowtie_r2,
        bt2idx="{query}/bowtie/{query}_chunk{chunk_num}.1.bt2l"
    log:
        "{query}/bam/{sample}_chunk{chunk_num}.log"
    output:
        bam_file=temp("{query}/bam/PE_{sample}_sorted_chunk{chunk_num}.bam")
    benchmark:
        repeat("benchmarks/bowtie_alignment_{query}_{sample}_chunk{chunk_num}.benchmark.txt", 1)
    params:
        index="{query}/bowtie/{query}_chunk{chunk_num}"
    threads:
        cpu_count()
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -1 {input.fastq_r1} -2 {input.fastq_r2} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"



def get_sorted_bam_paths(wildcards):

    get_paths = checkpoints.count_bt2_idx.get(query=wildcards.query)
    idx_chunk_total = len(pd.read_csv(get_paths.output[0], sep='\t'))

    if PE_MODERN:
        return expand("{query}/bam/{reads}_{sample}_sorted_chunk{chunk_num}.bam", query=wildcards.query, reads=['PE'],
            sample=wildcards.sample, chunk_num=[x+1 if idx_chunk_total > 1 else 1 for x in range(idx_chunk_total)] )
    if PE_ANCIENT or SE:
        return expand("{query}/bam/{reads}_{sample}_sorted_chunk{chunk_num}.bam", query=wildcards.query, reads=['SE'],
            sample=wildcards.sample, chunk_num=[x+1 if idx_chunk_total > 1 else 1 for x in range(idx_chunk_total)])



rule merge_bams_single_end:
    input:
        # aln_done="{query}/bam/{sample}_all_aln.done",
        aln_path=get_sorted_bam_paths
    log:
        "{query}/bam/{sample}_merge_bams.log"
    output:
        "{query}/bam/SE_{sample}_sorted.bam"
    shell:
        "samtools merge -f {output} {input.aln_path} 2> {log}"



rule merge_bams_paired_end:
    input:
        # aln_done="{query}/bam/{sample}_all_aln.done",
        aln_path=get_sorted_bam_paths
    log:
        "{query}/bam/{sample}_merge_bams.log"
    output:
        "{query}/bam/PE_{sample}_sorted.bam"
    shell:
        "samtools merge -f {output} {input.aln_path} 2> {log}"



# rule dedup_merged:
#     input:
#         "{query}/bam/{sample}_sorted.bam"
#     log:
#         "{query}/bam/{sample}_sorted_rmdup.log"
#     output:
#         "{query}/bam/{sample}_sorted_rmdup.bam"
#     benchmark:
#         repeat("benchmarks/dedup_merged_{query}_{sample}.benchmark.txt", 1)
#     params:
#         output="{query}/bam/"
#     shell:
#         "dedup --merged --input {input} --output {params.output} &> {log}"



rule extract_fastq_single_end:
    input:
        "{query}/bam/SE_{sample}_sorted.bam"
    log:
        "{query}/fastq/SE/{sample}_mapq.log"
    output:
        "{query}/fastq/SE/{sample}_mapq.fastq.gz"
    benchmark:
        repeat("benchmarks/extract_fastq_single_end_{query}_{sample}.benchmark.txt", 3)
    params:
        min_mapq=config['min_mapq']
    shell:
        # "( samtools view -h -F 4 {input} | samtools fastq -c 6 - > {output} ) 2> {log}"
        "( samtools view -h -F 4 {input} | samtools fastq -c 6 - | seqkit rmdup -n -o {output} ) 2> {log}"



rule extract_fastq_paired_end:
    input:
        "{query}/bam/PE_{sample}_sorted.bam"
    log:
        "{query}/fastq/PE/{sample}_mapq.log"
    output:
        "{query}/fastq/PE/{sample}_R1_mapq.fastq.gz",
        "{query}/fastq/PE/{sample}_R2_mapq.fastq.gz"
    benchmark:
        repeat("benchmarks/extract_fastq_paired_end_{query}_{sample}.benchmark.txt", 3)
    params:
        min_mapq=config['min_mapq']
    shell:
        # "( samtools view -h -F 4 {input} "
        # "| samtools fastq -c 6 -1 {output[0]} -2 {output[1]} -0 /dev/null -s /dev/null - ) 2> {log}"
        "( samtools view -h -F 4 {input} "
        "| samtools fastq -c 6 -1 {wildcards.query}/fastq/PE/temp_R1.fastq.gz "
        "-2 {wildcards.query}/fastq/PE/temp_R2.fastq.gz -0 /dev/null -s /dev/null - );"
        "seqkit rmdup -n {wildcards.query}/fastq/PE/temp_R1.fastq.gz -o {output[0]}; "
        "seqkit rmdup -n {wildcards.query}/fastq/PE/temp_R2.fastq.gz -o {output[1]};"
        "rm {wildcards.query}/fastq/PE/temp_R1.fastq.gz; rm {wildcards.query}/fastq/PE/temp_R2.fastq.gz  2> {log}"


# ruleorder: extract_fastq_paired_end > extract_fastq_single_end



rule average_fastq_read_len_single_end:
    input:
        "{query}/fastq/SE/{sample}_mapq.fastq.gz"
    log:
        "{query}/fastq/SE/{sample}_mapq_readlen.log"
    output:
        "{query}/fastq/SE/{sample}_mapq.readlen"
    benchmark:
        repeat("benchmarks/average_fastq_read_len_single_end_{query}_{sample}.benchmark.txt", 3)
    params:
        sample_size=SUBSAMPLE_FIXED_READS
    shell:
        "seqtk sample {input} {params.sample_size} | seqtk seq -A | grep -v '^>' | "
        "awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output} 2> {log}"



rule average_fastq_read_len_paired_end:
    input:
        mate1="{query}/fastq/PE/{sample}_R1_mapq.fastq.gz",
        mate2="{query}/fastq/PE/{sample}_R2_mapq.fastq.gz"
    log:
        "{query}/fastq/PE/{sample}_mapq_readlen.log"
    output:
        mate1=temp("{query}/fastq/{sample}_R1_mapq.readlen"),
        mate2=temp("{query}/fastq/{sample}_R2_mapq.readlen"),
        pair="{query}/fastq/PE/{sample}_mapq_pair.readlen"
    benchmark:
        repeat("benchmarks/average_fastq_read_len_paired_end_{query}_{sample}.benchmark.txt", 3)
    params:
        sample_size=SUBSAMPLE_FIXED_READS
    shell:
        "seqtk sample {input.mate1} {params.sample_size} | seqtk seq -A | grep -v '^>' | "
        "awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output.mate1} 2> {log}; "
        "seqtk sample {input.mate2} {params.sample_size} | seqtk seq -A | grep -v '^>' | "
        "awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output.mate2} 2> {log}; "
        "cat {output.mate1} {output.mate2} | awk '{{sum += $1; n++ }} END {{if (n>0) print sum/n;}}' "
        "1> {output.pair} 2> {log} "
