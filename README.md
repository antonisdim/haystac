RIP 
===

Introduction 
------------

Rip is a pipeline for species identification from singe (sample DNA comes from one organism) or multiple (sample DNA is a metagenome) source samples.

Quick Start 
-----------

Four commands/modules that allow you to run the method from start to finish 

```
conda install -c bioconda rip
rip_multilevel database --query '("Yersinia"[Organism] OR "Yersinia"[Organism]) AND "complete genome"[All Fields]' --db-output ./yersinia_example
rip_multilevel sample --sra SRR12157896 --collapse --sample-output-dir samples
rip_multilevel analyse --mode abundances --database yersinia_example/database_build_config.yaml --sample ./samples/SRR12157896_config.yaml --analysis-output-dir ./analysis_output
```

Alignment based metagenomics
----------------------------

Rip is an easy to use pipeine for metagenomic identifications 

You can easily:

1. Construct a database 
2. Prepare your sample for analysis, including trimming sequencing adapters, collapsing paired end reads and also downloading data from the SRA
3. Perform species identification 
4. Perform a chemical damage pattern analysis 

Setup 
-----

Rip can be run on any unix based system. It needs the conda package manager to run and it easy to install. 

```example commmands when I got them```


License 
-------

MIT

Citation 
--------

Happening soon 
