#!/usr/bin/env python
# -*- coding: utf-8 -*-

from multiprocessing import cpu_count

##### Target rules #####


rule bowtie_index:
    input:
         "{query}/bowtie/{query}.fasta"
    log:
         "{query}/bowtie/{query}.index.log"
    output:
         expand("{{query}}/bowtie/{{query}}.{n}.bt2", n=[1, 2, 3, 4]),
         expand("{{query}}/bowtie/{{query}}.rev.{n}.bt2", n=[1, 2])
    shell:
          "bowtie2-build {input} {wildcards.query}/bowtie/{wildcards.query} &> {log}"


rule bowtie_alignment:
    input:
        fastq=lambda wildcards: config['samples'][wildcards.sample],
        bt2idx="{query}/bowtie/{query}.1.bt2"
    log:
        "{query}/bam/{sample}.log"
    params:
        index="{query}/bowtie/{query}",
    output:
        "{query}/bam/{sample}_sorted.bam"
    threads:
        cpu_count()
    shell:
         "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -U {input.fastq} "
         "| samtools sort -O bam -o {output} ) 2> {log}"


rule remove_duplicates:
    input:
        "{query}/bam/{sample}_sorted.bam"
    log:
        "{query}/bam/{sample}_sorted_rmdup.log"
    output:
        "{query}/bam/{sample}_sorted_rmdup.bam"
    params:
        output="{query}/bam/"
    shell:
        "dedup --merged --input {input} --output {params.output} &> {log}"


rule extract_fastq:
    input:
        "{query}/bam/{sample}_sorted_rmdup.bam"
    log:
        "{query}/bam/{sample}_fastq.log"
    output:
        "{query}/bam/{sample}.fastq"
    params:
        config['alignment_qscore']
    shell:
        "samtools view -h -q {params} {input} | samtools fastq - > {output} &> {log}"


rule count_fastq_length:
    input:
         fastq=lambda wildcards: config['samples'][wildcards.sample]
    log:
         "fastq/{sample}.log"
    output:
         "fastq/{sample}.size"
    shell:
          "expr $(gunzip -c {input.fastq} | wc -l) / 4 1> {output} 2> {log}"

