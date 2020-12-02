HAYSTAC: High-AccuracY and Scalable Taxonomic Assignment of MetagenomiC data 
===

Introduction 
------------

Alignment based metagenomics
----------------------------

`haystack` is an easy to use pipeline for metagenomic identifications 

You can easily:

1. Construct a database 
2. Prepare your sample for analysis, including trimming sequencing adapters, collapsing paired end reads and also downloading data from the SRA
3. Perform species identification 
4. Perform a chemical damage pattern analysis 

Setup 
-----

`haystac` can be run on any unix based system. It needs the conda package manager to run and it easy to install. 

Quick Start 
-----------

Four commands/modules that allow you to run the method from start to finish 

```
# install haystack from conda
conda install -c bioconda haystac
```

```
# build a target database of species you are interested in
haystac database \
    --query '("Yersinia"[Organism] OR "Yersinia"[Organism]) AND "complete genome"[All Fields]' \
    --output yersinia_db
```

```
# prepare a sample for analyis 
haystac sample \
    --sra SRR12157896 \
    --collapse \
    --output SRR12157896
```

```
# analyse a sample against a specific database
haystac analyse \
    --mode abundances \
    --database yersinia_db \
    --sample SRR12157896 \
    --output yersinia_SRR12157896
```

```example commmands when I got them```


License 
-------

MIT

Citation 
--------

Happening soon 
