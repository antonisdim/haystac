#!/usr/bin/env python
# -*- coding: utf-8 -*-

from multiprocessing import cpu_count

SUBSAMPLE_FIXED_READS = 200000
WITH_REFSEQ_REP = config['WITH_REFSEQ_REP']
SRA_LOOKUP = config['SRA_LOOKUP']
PE_ANCIENT = config['PE_ANCIENT']
PE_MODERN = config['PE_MODERN']
SE = config['SE']


##### Target rules #####


def filtering_bowtie_aln_inputs(wildcards):
    if WITH_REFSEQ_REP:
        return ["{query}/bowtie/{query}_entrez.fasta.gz".format(query=wildcards.query),
                "{query}/bowtie/{query}_refseq_prok.fasta.gz".format(query=wildcards.query)]
    else:
        return ["{query}/bowtie/{query}_entrez.fasta.gz".format(query=wildcards.query)]



rule bowtie_index:
    input:
        filtering_bowtie_aln_inputs
    log:
        "{query}/bowtie/{query}_index.log"
    output:
        expand("{{query}}/bowtie/{{query}}.{n}.bt2l", n=[1, 2, 3, 4]),
        expand("{{query}}/bowtie/{{query}}.rev.{n}.bt2l", n=[1, 2])
    benchmark:
        repeat("benchmarks/bowtie_index_{query}.benchmark.txt", 1)
    run:
        if WITH_REFSEQ_REP:
            shell("cat {input} > {wildcards.query}/bowtie/{wildcards.query}.fasta.gz; "
                  "bowtie2-build --large-index {wildcards.query}/bowtie/{wildcards.query}.fasta.gz "
                  "{wildcards.query}/bowtie/{wildcards.query} &> {log}")
        else:
            shell("bowtie2-build --large-index {input} {wildcards.query}/bowtie/{wildcards.query} &> {log}")



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
            return config['samples'][wildcards.sample]['R1']
        elif PE_ANCIENT:
            return config['samples'][wildcards.sample]
        elif SE:
            return config['samples'][wildcards.sample]



def get_inputs_for_bowtie_r2(wildcards):
    if SRA_LOOKUP:
        if PE_MODERN:
            return "fastq_inputs/PE_mod/{sample}_R2_adRm.fastq.gz".format(sample=wildcards.sample)
        else:
            return ""

    else:
        if PE_MODERN:
            return config['samples'][wildcards.sample]['R2']
        else:
            return ""



rule bowtie_alignment:
    input:
        fastq=get_inputs_for_bowtie_r1,
        bt2idx="{query}/bowtie/{query}.1.bt2l"
    log:
        "{query}/bam/{sample}.log"
    params:
        index="{query}/bowtie/{query}",
    output:
        bam_file="{query}/bam/SE_{sample}_sorted.bam"
    benchmark:
        repeat("benchmarks/bowtie_alignment_{query}_{sample}.benchmark.txt", 1)
    threads:
        cpu_count()
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -U {input.fastq} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"



rule bowtie_alignment_paired_end:
    input:
        fastq_r1=get_inputs_for_bowtie_r1,
        fastq_r2=get_inputs_for_bowtie_r2,
        bt2idx="{query}/bowtie/{query}.1.bt2l"
    log:
        "{query}/bam/{sample}.log"
    output:
        bam_file="{query}/bam/PE_{sample}_sorted.bam"
    benchmark:
        repeat("benchmarks/bowtie_alignment_paired_end_{query}_{sample}.benchmark.txt", 1)
    params:
        index="{query}/bowtie/{query}"
    threads:
        cpu_count()
    shell:
        "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -1 {input.fastq_r1} -2 {input.fastq_r2} "
        "| samtools sort -O bam -o {output.bam_file} ) 2> {log}"



# todo this needs to be implemented differently in the pipeline
rule dedup_merged:
    input:
        "{query}/bam/{sample}_sorted.bam"
    log:
        "{query}/bam/{sample}_sorted_rmdup.log"
    output:
        "{query}/bam/{sample}_sorted_rmdup.bam"
    benchmark:
        repeat("benchmarks/dedup_merged_{query}_{sample}.benchmark.txt", 1)
    params:
        output="{query}/bam/"
    shell:
        "dedup --merged --input {input} --output {params.output} &> {log}"



rule extract_fastq_single_end:
    input:
        "{query}/bam/SE_{sample}_sorted.bam"
    log:
        "{query}/fastq/{sample}_mapq.log"
    output:
        "{query}/fastq/SE/{sample}_mapq.fastq.gz"
    benchmark:
        repeat("benchmarks/extract_fastq_single_end_{query}_{sample}.benchmark.txt", 3)
    params:
        min_mapq=config['min_mapq']
    shell:
        "( samtools view -h -F 4 {input} | samtools fastq -c 6 - > {output} ) 2> {log}"



rule extract_fastq_paired_end:
    input:
        "{query}/bam/PE_{sample}_sorted.bam"
    log:
        "{query}/fastq/{sample}_mapq.log"
    output:
        "{query}/fastq/PE/{sample}_R1_mapq.fastq.gz",
        "{query}/fastq/PE/{sample}_R2_mapq.fastq.gz"
    benchmark:
        repeat("benchmarks/extract_fastq_paired_end_{query}_{sample}.benchmark.txt", 3)
    params:
        min_mapq=config['min_mapq']
    shell:
        "( samtools view -h -F 4 {input} "
        "| samtools fastq -c 6 -1 {output[0]} -2 {output[1]} -0 /dev/null -s /dev/null - ) 2> {log}"



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
