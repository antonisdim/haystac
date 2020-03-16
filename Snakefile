#!/usr/bin/env python
# -*- coding: utf-8 -*-

configfile: "config.yaml"

##### Modules #####

include: "rules/entrez.smk"
include: "rules/bowtie.smk"
include: "rules/sigma.smk"  # TODO rename file now that we're not using sigma
include: "rules/metagenomics.smk"

##### Target rules #####

rule all:
    input:
        # "yersinia_test/entrez/yersinia_test-nuccore.tsv",
        # "yersinia_test/entrez/yersinia_test-taxa.tsv",
        # "yersinia_test/bowtie/yersinia_test.fasta",
        # "yersinia_test/bowtie/yersinia_test.1.bt2l",
        # # "yersinia_test/bam/RISE00_sorted.bam",
        # "yersinia_test/bam/RISE00_sorted_rmdup.bam",
        # "yersinia_test/fastq/RISE00_mapq.fastq",
        # "fastq/RISE00.size",
        # "yersinia_test/fastq/RISE00_mapq.readlen",
        # expand("database/{taxon}/{accession}_index.done", zip,
        #         taxon=["Yersinia_aldovae", "Yersinia_aleksiciae", "Yersinia_canariae", "Yersinia_enterocolitica", "Yersinia_frederiksenii", "Yersinia_hibernica", "Yersinia_intermedia", "Yersinia_kristensenii", "Yersinia_massiliensis", "Yersinia_pestis", "Yersinia_pseudotuberculosis", "Yersinia_rohdei", "Yersinia_ruckeri", "Yersinia_similis", "Yersinia_sp._FDAARGOS_228", "Yersinia_sp._KBS0713"],
        #         accession=["NZ_CP009781.1", "NZ_CP011975.1", "NZ_CP043727.1", "NZ_HF571988.1", "NZ_CP023962.1", "NZ_CP032487.1", "NZ_CP027397.1", "NZ_CP009997.1", "NZ_CP028487.1", "NZ_CP033699.1", "NZ_CP033713.1", "NZ_CP009787.1", "NZ_CP017236.1", "NZ_CP007230.1", "NZ_CP020409.2", "NZ_CP042173.1"]),
        # expand("yersinia_test/sigma/RISE00/{taxon}/{taxon}_{accession}.bam", zip,
        #        taxon=["Yersinia_aldovae", "Yersinia_aleksiciae", "Yersinia_canariae", "Yersinia_enterocolitica", "Yersinia_frederiksenii", "Yersinia_hibernica", "Yersinia_intermedia", "Yersinia_kristensenii", "Yersinia_massiliensis", "Yersinia_pestis", "Yersinia_pseudotuberculosis", "Yersinia_rohdei", "Yersinia_ruckeri", "Yersinia_similis", "Yersinia_sp._FDAARGOS_228", "Yersinia_sp._KBS0713"],
        #        accession=["NZ_CP009781.1", "NZ_CP011975.1", "NZ_CP043727.1", "NZ_HF571988.1", "NZ_CP023962.1", "NZ_CP032487.1", "NZ_CP027397.1", "NZ_CP009997.1", "NZ_CP028487.1", "NZ_CP033699.1", "NZ_CP033713.1", "NZ_CP009787.1", "NZ_CP017236.1", "NZ_CP007230.1", "NZ_CP020409.2", "NZ_CP042173.1"]),
        # expand("yersinia_test/ts_tv_counts/RISE00/{taxon}_count_{accession}.csv", zip,
        #        taxon=["Yersinia_aldovae", "Yersinia_aleksiciae", "Yersinia_canariae", "Yersinia_enterocolitica", "Yersinia_frederiksenii", "Yersinia_hibernica", "Yersinia_intermedia", "Yersinia_kristensenii", "Yersinia_massiliensis", "Yersinia_pestis", "Yersinia_pseudotuberculosis", "Yersinia_rohdei", "Yersinia_ruckeri", "Yersinia_similis", "Yersinia_sp._FDAARGOS_228", "Yersinia_sp._KBS0713"],
        #        accession=["NZ_CP009781.1", "NZ_CP011975.1", "NZ_CP043727.1", "NZ_HF571988.1", "NZ_CP023962.1", "NZ_CP032487.1", "NZ_CP027397.1", "NZ_CP009997.1", "NZ_CP028487.1", "NZ_CP033699.1", "NZ_CP033713.1", "NZ_CP009787.1", "NZ_CP017236.1", "NZ_CP007230.1", "NZ_CP020409.2", "NZ_CP042173.1"]),
        # expand("yersinia_test/ts_tv_counts/RISE00/{taxon}_count_{accession}_aggregating.done", zip,
        #        taxon=["Yersinia_aldovae", "Yersinia_aleksiciae", "Yersinia_canariae", "Yersinia_enterocolitica", "Yersinia_frederiksenii", "Yersinia_hibernica", "Yersinia_intermedia", "Yersinia_kristensenii", "Yersinia_massiliensis", "Yersinia_pestis", "Yersinia_pseudotuberculosis", "Yersinia_rohdei", "Yersinia_ruckeri", "Yersinia_similis", "Yersinia_sp._FDAARGOS_228", "Yersinia_sp._KBS0713"],
        #        accession=["NZ_CP009781.1", "NZ_CP011975.1", "NZ_CP043727.1", "NZ_HF571988.1", "NZ_CP023962.1", "NZ_CP032487.1", "NZ_CP027397.1", "NZ_CP009997.1", "NZ_CP028487.1", "NZ_CP033699.1", "NZ_CP033713.1", "NZ_CP009787.1", "NZ_CP017236.1", "NZ_CP007230.1", "NZ_CP020409.2", "NZ_CP042173.1"]),
        # "yersinia_test/probabilities/RISE00/TsTvMatrix.csv",
        # "yersinia_test/probabilities/RISE00/Posterior_Probabilities.csv",
        #  expand("yersinia_test/probabilities/RISE00/{taxon}_t_test_pvalue_{accession}.txt", zip,
        #        taxon=["Yersinia_aldovae", "Yersinia_aleksiciae", "Yersinia_canariae", "Yersinia_enterocolitica", "Yersinia_frederiksenii", "Yersinia_hibernica", "Yersinia_intermedia", "Yersinia_kristensenii", "Yersinia_massiliensis", "Yersinia_pestis", "Yersinia_pseudotuberculosis", "Yersinia_rohdei", "Yersinia_ruckeri", "Yersinia_similis", "Yersinia_sp._FDAARGOS_228", "Yersinia_sp._KBS0713"],
        #        accession=["NZ_CP009781.1", "NZ_CP011975.1", "NZ_CP043727.1", "NZ_HF571988.1", "NZ_CP023962.1", "NZ_CP032487.1", "NZ_CP027397.1", "NZ_CP009997.1", "NZ_CP028487.1", "NZ_CP033699.1", "NZ_CP033713.1", "NZ_CP009787.1", "NZ_CP017236.1", "NZ_CP007230.1", "NZ_CP020409.2", "NZ_CP042173.1"]),
        #  "yersinia_test/probabilities/RISE00/t_test_pvalues.txt",
        #  "yersinia_test/probabilities/RISE00/RISE00_posterior_abundance.tsv"
        "example1/entrez/example1-nuccore.tsv",
        "example1/entrez/example1-taxa.tsv",
        "example1/bowtie/example1.fasta.gz"




