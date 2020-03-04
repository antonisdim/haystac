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
         expand("{{query}}/bowtie/{{query}}.{n}.bt2l", n=[1, 2, 3, 4]),
         expand("{{query}}/bowtie/{{query}}.rev.{n}.bt2l", n=[1, 2])
    shell:
          # TODO did you forget to tell bowtie2 to use a --large-index?
          "bowtie2-build {input} {wildcards.query}/bowtie/{wildcards.query} &> {log}"


rule bowtie_alignment:
    input:
        fastq=lambda wildcards: config['samples'][wildcards.sample],
        bt2idx="{query}/bowtie/{query}.1.bt2l"
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
        "{query}/fastq/{sample}_mapq.log"
    output:
        # TODO gzip the output
        "{query}/fastq/{sample}_mapq.fastq"
    shell:
        "( samtools view -h -q {config.min_mapq} {input} "
        "| samtools fastq - > {output} ) 2> {log}"


# TODO this rule has nothing to do with bowtie, so move it to another file where the output is actually used
rule count_fastq_length:
    input:
         fastq=lambda wildcards: config['samples'][wildcards.sample]
    log:
         "fastq/{sample}.log"
    output:
         "fastq/{sample}.size"
    shell:
          # TODO this should work for both gzipped and raw fastq input, maybe switch to using `seqtk`
          "expr $(gunzip -c {input.fastq} | wc -l) / 4 1> {output} 2> {log}"


rule average_fastq_read_len:
    input:
        "{query}/fastq/{sample}_mapq.fastq"
    log:
        "{query}/fastq/{sample}_mapq_readlen.log"
    output:
        "{query}/fastq/{sample}_mapq.readlen"
    shell:
         # TODO don't use magic numbers like 2000000, this should be a static constant
         # TODO use `seqtk sample` to subsample rather than head, as this fastq is sorted by species name
         #      e.g. `seqtk sample fastq/CSL.all.fastq.gz {size} | seqtk seq -A | grep -v '^>' | awk '{{count++; bases += length}} END{{print bases/count}}'`
         #      but you should also threshold this, as `seqtk sample` will be slow if size is much larger than the file itself
         """head -n 2000000 {input} | awk "{{if(NR%4==2) {{count++; bases += length}} }} END{{print bases/count}}" \
         1> {output} 2> {log}"""
