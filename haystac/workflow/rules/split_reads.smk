#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from haystac.workflow.scripts.utilities import get_total_paths, PE
from haystac.cli import CODE_DIR


rule get_dirichlet_reads:
    input:
        bam_file=config["analysis_output_dir"] + "/alignments/{sample}/{reads}/{orgname}/{orgname}_{accession}.bam",
        dirichlet_matrix=(
            config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv"
        ),
    output:
        config[
            "analysis_output_dir"
        ] + "/dirichlet_reads/{sample}/{orgname}/{orgname}_{accession}_dirichlet_{reads}.bam",
    benchmark:
        repeat(
            "benchmarks/get_dirichlet_reads_{sample}_{reads}_{orgname}_{accession}.benchmark.txt", 1,
        )
    message:
        "Preparing bam files with the Dirichlet assigned reads for taxon {wildcards.orgname} "
        "for sample {wildcards.sample} "
    script:
        "../scripts/get_dirichlet_reads.py"


rule get_grey_matter_reads_se:
    input:
        fastq=config["analysis_output_dir"] + "/fastq/SE/{sample}_mapq.fastq.gz",
        dirichlet_matrix=(
            config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv"
        ),
    output:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet.fastq.gz",
    benchmark:
        repeat("benchmarks/get_{sample}_Grey_Matter.benchmark.txt", 1)
    message:
        "Preparing fastq files with all the reads that got assigned to the Grey Matter for sample {wildcards.sample}."
    shell:
        "python {CODE_DIR}/workflow/scripts/get_matter_reads.py --input_fastq {input.fastq} "
        "--matrix_file {input.dirichlet_matrix} --output_fastq {output} --matter grey"


rule get_grey_matter_reads_pe:
    input:
        fastq_r1=config["analysis_output_dir"] + "/fastq/PE/{sample}_R1_mapq.fastq.gz",
        fastq_r2=config["analysis_output_dir"] + "/fastq/PE/{sample}_R2_mapq.fastq.gz",
        dirichlet_matrix=(
            config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv"
        ),
    output:
        out_r1=(
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_R1.fastq.gz"
        ),
        out_r2=(
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_R2.fastq.gz"
        ),
    benchmark:
        repeat("benchmarks/get_{sample}_Grey_Matter.benchmark.txt", 1)
    message:
        "Preparing fastq files with all the reads that got assigned to the Grey Matter for sample {wildcards.sample}."
    shell:
        "python {CODE_DIR}/workflow/scripts/get_matter_reads.py --input_fastq {input.fastq_r1} "
        "--matrix_file {input.dirichlet_matrix} "
        "--output_fastq {output.out_r1} --matter grey; python ../scripts/get_matter_reads.py "
        "--input_fastq {input.fastq_r2} --matrix_file {input.dirichlet_matrix} --output_fastq {output.out_r2} "
        "--matter grey"


rule get_dark_matter_reads_se:
    input:
        fastq=config["fastq"] if config["fastq"] else config["fastq_r1"],
        dirichlet_matrix=(
            config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv"
        ),
    output:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet.fastq.gz",
    benchmark:
        repeat("benchmarks/get_{sample}_Dark_Matter.benchmark.txt", 1)
    message:
        "Preparing fastq files with all the reads that got assigned to the Dark Matter for sample {wildcards.sample}."
    shell:
        "python {CODE_DIR}/workflow/scripts/get_matter_reads.py --input_fastq {input.fastq} "
        "--matrix_file {input.dirichlet_matrix} "
        "--output_fastq {output} --matter dark"


rule get_dark_matter_reads_pe:
    input:
        fastq_r1=config["fastq_r1"] if config["fastq_r1"] else config["fastq"],
        fastq_r2=config["fastq_r2"],
        dirichlet_matrix=(
            config["analysis_output_dir"] + "/probabilities/{sample}/{sample}_likelihood_ts_tv_matrix.csv"
        ),
    output:
        out_r1=(
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet_R1.fastq.gz"
        ),
        out_r2=(
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet_R2.fastq.gz"
        ),
    benchmark:
        repeat("benchmarks/get_{sample}_Dark_Matter.benchmark.txt", 1)
    message:
        "Preparing fastq files with all the reads that got assigned to the Dark Matter for sample {wildcards.sample}."
    shell:
        "python {CODE_DIR}/workflow/scripts/get_matter_reads.py --input_fastq {input.fastq_r1} "
        "--matrix_file {input.dirichlet_matrix} "
        "--output_fastq {output.out_r1} --matter dark; python ../scripts/get_matter_reads.py "
        "--input_fastq {input.fastq_r2} --matrix_file {input.dirichlet_matrix} --output_fastq {output.out_r2} "
        "--matter dark"


def get_dirichlet_bams(_):
    """Get paths for dirichlet assigned reads in bams"""
    return [
        config["analysis_output_dir"]
        + "/dirichlet_reads/{sample}/"
        + f"{orgname}/{orgname}_{accession}_dirichlet_"
        + config["read_mode"]
        + ".bam"
        for orgname, accession in get_total_paths(checkpoints, config)
    ]


rule all_dirichlet:
    input:
        get_dirichlet_bams,
        [
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet_R1.fastq.gz"
            if config["read_mode"] == PE
            else config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Grey_Matter/Grey_Matter_dirichlet.fastq.gz"
        ],
        [
            config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet_R1.fastq.gz"
            if config["read_mode"] == PE
            else config["analysis_output_dir"] + "/dirichlet_reads/{sample}/Dark_Matter/Dark_Matter_dirichlet.fastq.gz"
        ],
    output:
        config["analysis_output_dir"] + "/dirichlet_reads/{sample}_dirichlet_reads.done",
    benchmark:
        repeat("benchmarks/all_dirichlet_reads_{sample}.benchmark.txt", 1)
    message:
        "All the dirichlet assigned reads have been put in bam files for {wildcards.sample}."
    shell:
        "touch {output}"
