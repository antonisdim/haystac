#!/usr/bin/env python
# -*- coding: utf-8 -*-

from multiprocessing import cpu_count

SUBSAMPLE_FIXED_READS = 200000

##### Target rules #####


rule bowtie_index:
    input:
         "{query}/bowtie/{query}.fasta.gz"
    log:
         "{query}/bowtie/{query}_index.log"
    output:
         expand("{{query}}/bowtie/{{query}}.{n}.bt2l", n=[1, 2, 3, 4]),
         expand("{{query}}/bowtie/{{query}}.rev.{n}.bt2l", n=[1, 2])
    benchmark:
        repeat("benchmarks/bowtie_index_{query}.benchmark.txt", 3)
    shell:
          "bowtie2-build --large-index {input} {wildcards.query}/bowtie/{wildcards.query} &> {log}"



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
    benchmark:
        repeat("benchmarks/bowtie_alignment_{query}_{sample}.benchmark.txt", 3)
    threads:
        cpu_count()
    shell:
         "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -U {input.fastq} "
         "| samtools sort -O bam -o {output} ) 2> {log}"



# todo it's provisional haven't used it
rule bowtie_alignment_paired_end:
    input:
        fastq_r1=lambda wildcards: config['samples'][wildcards.sample],
        fastq_r2=lambda wildcards: config['samples'][wildcards.sample],
        bt2idx="{query}/bowtie/{query}.1.bt2l",
    log:
        "{query}/bam/{sample}.log"
    output:
        bam_file="{query}/bam/{sample}_sorted.bam"
    benchmark:
        repeat("benchmarks/bowtie_alignment_paired_end_{query}_{sample}.benchmark.txt", 3)
    params:
        index="{query}/bowtie/{query}"
    threads:
        cpu_count()
    shell:
         "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -1 {input.fastq_r1} -2 {input.fastq_r2} "
         "| samtools sort -O bam -o {output} ) 2> {log}"



ruleorder: bowtie_alignment > bowtie_alignment_paired_end



rule dedup_merged:
    input:
        "{query}/bam/{sample}_sorted.bam"
    log:
        "{query}/bam/{sample}_sorted_rmdup.log"
    output:
        "{query}/bam/{sample}_sorted_rmdup.bam"
    benchmark:
        repeat("benchmarks/dedup_merged_{query}_{sample}.benchmark.txt", 3)
    params:
        output="{query}/bam/"
    shell:
        "dedup --merged --input {input} --output {params.output} &> {log}"



rule extract_fastq_single_end:
    input:
        "{query}/bam/{sample}_sorted_rmdup.bam"
    log:
        "{query}/fastq/{sample}_mapq.log"
    output:
        "{query}/fastq/{sample}_mapq.fastq.gz"
    benchmark:
        repeat("benchmarks/extract_fastq_single_end_{query}_{sample}.benchmark.txt", 3)
    params:
        min_mapq = config['min_mapq']
    shell:
        "( samtools view -h -q {params.min_mapq} {input} "
        "| samtools fastq -c 6 - > {output} ) 2> {log}"



rule extract_fastq_paired_end:
    input:
        "{query}/bam/{sample}_sorted_rmdup.bam"
    log:
        "{query}/fastq/{sample}_mapq.log"
    output:
        "{query}/fastq/{sample}_r1_mapq.fastq.gz",
        "{query}/fastq/{sample}_r2_mapq.fastq.gz"
    benchmark:
        repeat("benchmarks/extract_fastq_paired_end_{query}_{sample}.benchmark.txt", 3)
    params:
        min_mapq = config['min_mapq']
    shell:
        "( samtools view -h -q {params.min_mapq} {input} "
        "| samtools fastq -c 6 -1 {output[0]} -2 {output[1]} -0 /dev/null -s /dev/null - ) 2> {log}"



ruleorder: extract_fastq_single_end > extract_fastq_paired_end



rule average_fastq_read_len_single_end:
    input:
        "{query}/fastq/{sample}_mapq.fastq.gz"
    log:
        "{query}/fastq/{sample}_mapq_readlen.log"
    output:
        "{query}/fastq/{sample}_mapq.readlen"
    benchmark:
        repeat("benchmarks/average_fastq_read_len_single_end_{query}_{sample}.benchmark.txt", 3)
    params:
        sample_size=SUBSAMPLE_FIXED_READS
    shell:
         "seqtk sample {input} {params.sample_size} | seqtk seq -A | grep -v '^>' | "
         "awk '{{count++; bases += length}} END{{print bases/count}}' 1> {output} 2> {log}"



rule average_fastq_read_len_paired_end:
    input:
        mate1 = "{query}/fastq/{sample}_r1_mapq.fastq.gz",
        mate2 = "{query}/fastq/{sample}_r2_mapq.fastq.gz"
    log:
        "{query}/fastq/{sample}_mapq_readlen.log"
    output:
        mate1 = temp("{query}/fastq/{sample}_mapq_mate1.readlen"),
        mate2 = temp("{query}/fastq/{sample}_mapq_mate2.readlen"),
        pair = "{query}/fastq/{sample}_mapq.readlen"
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



ruleorder: average_fastq_read_len_single_end > average_fastq_read_len_paired_end