#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from haystac.cli import CODE_DIR
from haystac.workflow.scripts.utilities import get_final_db_paths, PE


rule get_dirichlet_reads:
    input:
        bam_file=config["analysis_output_dir"] + "/alignments/{sample}/{reads}/{orgname}/{orgname}_{accession}.bam",
        dirichlet_matrix=(
            config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv"
        ),
    log:
        config["analysis_output_dir"]
        +"/dirichlet_reads/{sample}/{orgname}/{orgname}_{accession}_dirichlet_{reads}.log",
    output:
        bam_out=(
            config["analysis_output_dir"]
            + "/dirichlet_reads/{sample}/{orgname}/{orgname}_{accession}_dirichlet_{reads}.bam"
        ),
        read_list=temp(
            config["analysis_output_dir"]
            + "/dirichlet_reads/{sample}/{orgname}/{orgname}_{accession}_dirichlet_{reads}_list.txt"
        ),
    message:
        "Preparing bam files with the Dirichlet assigned reads for taxon {wildcards.orgname} "
        "for sample {wildcards.sample}"
    conda:
        "../envs/picard.yaml"
    shell:
        "awk -F, '$1 == \"{wildcards.orgname}\"' {input.dirichlet_matrix} | "
        "awk -F, '$7 == \"1.0\" {{print $2}}' > {output.read_list};  "
        "if [ $(cat {output.read_list}  | wc -l) == '0' ]; "
        "then samtools view -b {input.bam_file} -H -o {output.bam_out}; "
        "else (picard FilterSamReads --INPUT {input.bam_file} --OUTPUT {output.bam_out} "
        "--READ_LIST_FILE {output.read_list} --SORT_ORDER coordinate --FILTER includeReadList) 2> {log}; "
        "fi"


def get_sample_bam(wildcards):
    """Function to get the correct bam file for matters"""

    return config["analysis_output_dir"] + "/bam/" + config["read_mode"] + f"_{wildcards.sample}_sorted.bam"


rule get_grey_matter_reads:
    input:
        bam=get_sample_bam,
        dirichlet_matrix=(
            config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv"
        ),
    log:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_list.log",
    output:
        read_list=temp(
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_list.txt"
        ),
        bam_out=config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet.bam",
    message:
        "Preparing fastq files with all the reads that got assigned to the Grey Matter for sample {wildcards.sample}."
    conda:
        "../envs/picard.yaml"
    shell:
        "awk -F, 'NR>1 {{arr[$2]+=$7}} END {{for (a in arr) print a, arr[a]}}' {input.dirichlet_matrix} | "
        "awk '$2 == 0 {{print $1}}' > {output.read_list};"
        "if [ $(cat {output.read_list}  | wc -l) == '0' ]; "
        "then samtools view -b {input.bam} -H -o {output.bam_out}; "
        "else (picard FilterSamReads --INPUT {input.bam} --OUTPUT {output.bam_out} "
        "--READ_LIST_FILE {output.read_list} --SORT_ORDER coordinate --FILTER includeReadList) 2> {log}; "
        "fi"


rule get_dark_matter_reads:
    input:
        bam=get_sample_bam,
    log:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet.log",
    output:
        temp_sample_bam=temp(
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/{sample}_pic_sort.bam"
        ),
        dark_bam=config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet.bam",
    message:
        "Preparing fastq files with all the reads that got assigned to the Dark Matter for sample {wildcards.sample}."
    conda:
        "../envs/picard.yaml"
    shell:
        "(picard SortSam --INPUT {input.bam} --OUTPUT {output.temp_sample_bam} --SORT_ORDER queryname; "
        "picard FilterSamReads --INPUT {output.temp_sample_bam} --OUTPUT {output.dark_bam} "
        "--SORT_ORDER coordinate --FILTER excludeAligned) 2> {log}"


def get_dirichlet_bams(wildcards):
    """Get paths for dirichlet assigned reads in bams"""
    return [
        config["analysis_output_dir"]
        + f"/dirichlet_reads/{wildcards.sample}/"
        + f"{orgname}/{orgname}_{accession}_dirichlet_"
        + config["read_mode"]
        + ".bam"
        for orgname, accession in get_final_db_paths(checkpoints)
    ]


rule all_dirichlet:
    input:
        get_dirichlet_bams,
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet.bam",
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet.bam",
    output:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}_dirichlet_reads.done",
    message:
        "All the dirichlet assigned reads have been put in bam files for {wildcards.sample}."
    shell:
        "touch {output}"
