HAYSTACK documentation
===============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Introduction
============

HAYSTAC is a comprehensive computational tool for identifying species from DNA sequence data. It can pre-process sequencing data (adapter trimming), build a database, analyse sequencing data and perform a deamination profile analysis. It can work in two different modes. One mode performs metagenomic identifications from samples containing multiple organisms, and outputs mean posterior species abundances. The second mode can perform species identifications from single organism samples and it outputs species assignment posterior probabilities.

Installation
============

HAYSTAC can be installed though pip, bioconda, and you can download the source from git as well.

Pip

    pip install haystac

Conda

    conda install -c bioconda haystac

Git 
    git clone haystac


Workflow
========

HAYSTAC was designed to be used as a species identifier for either single organism or metagenomic samples. 
The full execution of the pipeline includes the construction of a database, the processing of a sample file and lastly the analysis of a sample against a specific database. For that purpose three modules have been designed, each of which has its own outputs. 

First the `haystac database` module can be used to construct a database, based on the user's needs. The user must provide a path where the outputs of the database can be stored. 
Expected outputs for `haystac database`:
db_taxa_accessions.tsv: includes all the taxon/accession pairs that are included in a database. 
idx_database.done: file indicating that all the individual bowtie2 indices for each taxon have been prepared
entrez: directory containing all the results from querying the NCBI, including the nucleotide and taxonomy databases
bowtie: directory containing the bowtie2 index files for all out database, that will be used for the filtering alignment
database_inputs: directory containing the representative RefSeq table that is downloaded from NCBI. 

The `haystac sample` module prepares a sample to be analysed. It deals with trimming adapters and collapsing PE reads if needed and it counts the number of reads that are included in the sample file provided by the user. 
Expected outputs for `haystac sample`:
fastq_inputs: folder containing the outputs of the sample module.
	meta: directory that includes the read count file.
	SE/PE/COLLAPSED: directory containing the trimmed reads produced by AdapteRremoval

The `haystac analyse` module performs the analysis of a sample against a specific database. All its outputs are created under the path that the user specifies with the `--output` path. 
Expected outputs for `haystac analyse`:
bam: directory containing the bam file of the filtering alignment 
fastq: directory containing the filtered reads in fastq format along with their average read length 
alignments: directory where all the individual alignment bam files for each taxon in a database of the filtered reads are outputted
ts_tv_counts: directory where all the transition and transversion counts are stored per taxon
probabilities: directory where the likelihood matrix and final posterior abundances/probabilities are stored. The final output for abundance calculation has the suffix {sample-prefix}_posterior_abundance.tsv, 
mapdamage: directory that includes all the mapdamage profiles for every taxon in our Database
reads: directory including all the Dirichlet reads for each taxon in out database. 


Tutorial
========

Configuring HAYSTAC
---------------

If a user is running haystac for the first time on a machine, they might want to run the `haystac config` module first. The user can use this module to provide an API key for querying NCBI and their preferred path for the cache genomes folder, among other things. All of the options have default values that can be used. If a user later on wishes to change any of these parameters specifically they can either run `haystac config` to pass a value to a specific argument.

Here is an example command that allows the configuration of using conda as a package manager for running the other haystac modules.

    haystac config --use-conda True

Unless the user has a deep understanding of their dataset we advise to be cautious when changing the base mismatch probability that is used later on in the method's probabilistic model. We also advise caution when changing the bowtie2 file size scaling factor.

Building the database
---------------------

In order to build the database we will be using the `database` module of HAYSTAC.
First we need to know what organisms we would like to include in our database. Do we only need the complete genomes of a specific genus or do we want more genera? 

Constructing the Query
----------------------

After deciding what taxa we would like to include in our database, we need to construct an NCBI query that will return all the accessions that belong to the taxa that we are interested in. One of the best ways to construct such a query is to go on the website of  NCBI's Nucleotide database (https://www.ncbi.nlm.nih.gov/nucleotide/), type in our query and get the correctly formatted query string from the "Search details" box. We can then use that string for the construction of our database. 

For example if we would like to build a database of all the complete genomes of the species in the Yersinia genus we can use the following command:

    haystac database --mode build --query '("Yersinia"[Organism] OR "Yersinia"[Organism]) AND "complete genome"[All Fields]' --output yersinia_example

For each species (or any other user defined taxonomic rank), the longest sequence per taxon will be used to populate our database. 

Representative RefSeq species
-----------------------------

When constructing a database there is always the option to include the species of the representative RefSeq database as well. All you need to do is include the corresponding flag in your command. 

    haystack database --mode build --query '("Yersinia"[Organism] OR "Yersinia"[Organism]) AND "complete genome"[All Fields]' --output yersinia_example --refseq-rep

Providing custom accessions 
---------------------------

It is also possible to provide your own accessions for a selected species/taxon. For that you will need to prepare a tab delimited file with two columns. The first column is the name of the taxon, that cannot contain any special characters, other than an underscore ('_'), and the second column is a valid NCBI accession code. 

Here is an example of the contents of such a file:

    Yersinia_ruckeri	NZ_CP025800.1

The we can simply run the following command 

    haystack database --mode build --query '("Yersinia"[Organism] OR "Yersinia"[Organism]) AND "complete genome"[All Fields]' --output yersinia_example --accessions acc_example.txt

Providing custom sequences
--------------------------

It is also possible to provide your own sequences for a taxon. To do that you will need a a tab delimited file containing the the name of the taxon with no special characters, and an underscore ('_') instead of spaces, a user defined accession code and the path of the fasta file. The fasta file that the path point to can be either uncompressed or compressed with gzip/bgzip.

Here is an example of such a file:

    Yersinia_ruckeri	user_seq_1	~/example_sequence/user_seq_1.fasta

The we can simply run the following command 

    haystack database --mode build --query '("Yersinia"[Organism] OR "Yersinia"[Organism]) AND "complete genome"[All Fields]' --db-output yersinia_example --sequences seq_example.txt

Combinations
------------

All of the previous options can be combined into one command. It is important to note that only one sequence file per taxon is allowed in our database, and the priority goes user defined accessions or sequences, representative RefSeq and then user specified query. It should be noted that for custom NCBI queries, plasmids can be fetched if they are part of a genome assembly. The only exception to the one sequence file per taxon rule is plasmid sequences of the RefSeq representative that are not part of complete genome assemblies and for that reason they are downloaded separately. 

Database configuration yaml file
--------------------------------

After the database has been created a yaml file with all the configured and default options will be outputted in the ` --output` path that the user has provided. This file will be needed for the sample analysis as it contains all the database related information.

Index building 
--------------

For the first part of the analysis an index out of all the genomes that are included in our database needs to be build. This is a process that can take big amounts of memory depending on the number and the complexity of sequences that our database includes. For that reason the user can specify the desired amount of memory resources available to haystac and the program will try to build the required index. This can be specified through the `--mem` flag, that can be appended to the any of the commands shown above. Memory resources need to be specified in MB. If the memory resources provided are less than the size of the files that need to be indexed an error will be raised. 

Database building modes
-----------------------

For the complete construction of a database, sequences need to be downloaded and subsequently indexed.
By specifying `--mode build` to `haystac database`, the program downloads and indexes all the sequences that have been requested by the user in one step.
If a user would like to only download sequence data and index them later it is possible to do so, by specifying `haystac database --mode fetch`, to download the sequences first and then execute `haystac database --mode index` in order to perform the indexing.
If mode `fetch` is run first then mode `index` should be run subsequently, and not mode `build`, otherwise an error will be raised. 

Here is an example of building a database in two steps instead of one 

    haystac database --mode fetch --query '("Yersinia"[Organism] OR "Yersinia"[Organism]) AND "complete genome"[All Fields]' --output yersinia_example
    haystac database --mode index --output yersinia_example

Building a mitochondrial DNA database
-------------------------------------

When a user is providing a query about eukaryotes it is also possible to build a database with only mitochondrial genomes (by default whole genome assemblies will be fetched for a given query). In order to do that a user can specify the `--mtDNA` flag when running `haystack database`. We strongly advise against having a mixed database of full eukaryotic genome assemblies for certain taxa and only mtDNA sequences for other taxa, as this will bias the identifications towards the taxa with full genome assemblies.

Preparing a sample for analysis
-------------------------------

After our database is built we need to prepare our samples for analysis. For that purpose, we are using the `sample` module of haystac. The input files can be SE, PE or collapsed reads. If the reads are collapsed they are going to be treated as SE reads.

It is possible to trim sequencing adapters and collapse PE reads by specifying the relative flags. Samples (specific sequencing runs) can be also downloaded from the SRA if an sra run accession is provided. 

If you have SE or already collapsed reads you only need to specify a file path for the `--fastq` flag and a sample name thought the `--sample-prefix` flag. 
If your input is PE reads then you will need to specify file paths for both the `--fastq-r1` and `--fastq-r2` flags and a sample name thought the `--sample-prefix` flag. 
If you want to download files from the SRA all you need to do is provide an SRA accession for the `--sra` flag (the sample name will be the name of the SRA accession code therefore do not use that flag in that case).

Here is an example of downloading reads from the SRA, trimming sequencing adapters and collapsing reads. 

    haystac sample --sra ERR1018966 --output sample_example

Sample configuration yaml file
--------------------------------

After the sample preparation has finished a yaml file with all the configured and default options will be created in the ` --output` path that the user has provided. This file will be needed for the sample analysis, as it contains the necessary information for each sample.

Sample analysis
---------------

In order to analyse any sample we will need to use the `analyse` module of haystac. The yaml files that specify the database and sample related options will also be needed for the `--database` and`--sample` flags.

Filtering Alignment
-------------------

The first step for the sample analysis is to filter in all the reads that align to any of the genomes in our database. For that we will need to use the `haystac analyse --mode filter`.

Here is an example command 

    haystac analyse --mode filter --database yersinia_example --sample sample_example --output analysis_output

Database Alignments
-------------------

After we have filtered our libraries we can align the filtered reads against all the genomes that are included in our database. This can be done by using mode `align` of `haystac analyse`.

For example

    haystac analyse --mode align --database yersinia_example --sample sample_example --output analysis_output

Likelihood calculation
----------------------

After all the individual alignments have been competed, the number of transitions and transversions will be counted for every read that has aligned against any of the reference genomes in our database. Then the likelihoods and posterior probabilities for each read being sampled from a given reference genome will be calculated. For this step we can use the `likelihoods` mode of `haystac analyse`.

    haystack analyse --mode likelihoods --database yersinia_example --sample sample_example --output analysis_output

Important Note on the Dirichlet Assignment process during Likelihood calculation
--------------------------------------------------------------------------------

It is important to be aware of the individual read posterior probability threshold, for a read to be assigned to a taxon. As a default HAYSTACK uses the conservative 0.75 probability threshold for the Dirichlet assignment. The higher value you pick the more conservative the assignments become. It is useful to sometimes pick a value depending on what taxa are being identified. If there is a need to distinguish between closely related taxa then a more conservative threshold would increase the specificity of the analysis therefore being more appropriate, whereas when you're trying to generally characterise a metagenome a less conservative value could increase the sensitivity of the analysis be more helpful.

Single organism sample or metagenome ? 
--------------------------------------

Depending on whether we would like to identify the species a sample is belongs to, or perform a metagenomic analysis, we can use the `probabilities` or `abundances` mode of `haystack analyse` respectively.

Assignment Probability Calculation
----------------------------------

In order to calculate posterior assignment probabilities we can run the following command

    haystac analyse --mode probabilities --database yersinia_example --sample sample_example --output analysis_output

Mean Posterior Abundances
-------------------------

In order to calculate mean posterior abundances we can run the following command 

    haystac analyse --mode abundances --database yersinia_example --sample sample_example --output analysis_output
    
Reads
-----

After the mean posterior abundances have been calculated for a sample, all the reads that have been assigned to a taxon through the Dirichlet process can be outputted in separate bam files ready for further downstream analyses (like assembling or variant calling for instance) via the `haystac reads` module. Reads that have been assigned to the Grey and Dark Matter are outputted in fastq files as they have not been uniquely assigned to a taxon.

Here is an example command

    haystac analyse --mode reads --database yersinia_example --sample sample_example --output analysis_output

Mapdamage analysis
------------------

If our samples are ancient we can use mapDamage to estimate the level of deamination in the reads that have aligned to any taxon in our database. For that we can use the `mapdamage` module of haystack. The mapDamage analysis will be performed on the subset of reads that have been uniquely assigned to a taxon through the dirichlet process. This module can be either run independently or after the `haystack reads` module.

Here is an example command 

    haystac analyse --mode mapdamage --database yersinia_example --sample sample_example --output analysis_output

Important note on sample analysis
---------------------------------

The first 3 steps (filter, align, likelihoods) can be executed automatically when you call the `probabilities` or `abundances` mode of haystac.


Command Line Interface
======================

haystac config
--------------

    -h, --help            show this help message and exit
    -e , --email          Email address for NCBI identification. Mandatory.
    -gc , --cache
                          Path where all the genomes that are downloaded and/or
                          used by haystack are being stored. (default ~/haystack/cache/)
    -b , --batchsize      Batchsize for fetching records from NCBI <int>
                          (default: 5)
    -mp , --mismatch-probability 
                          Base mismatch probability <float> (default: 0.05)
    -s , --bowtie2-scaling 
                          Factor to rescale/chunk the input file for the
                          mutlifasta index for the filtering alignment (default:
                          2.5)

haystac database
----------------

    -h, --help            show this help message and exit
    --dry-run
    -m , --mode           Database creation mode for haystack
    -o , --db-output      Path to the database output directory.
    -R , --refseq-rep     Use the prokaryotic representative species of the
                          RefSeq DB for the species id pipeline. only species no
                          strains. either or both of with_refseq_rep and
                          with_entrez_query should be set (default: False)
    -MT , --mtDNA         Download mitochondrial genomes for eukaryotes only. Do
                          not use with --refseq-rep or any queries for
                          prokaryotes (default: False)
    -q , --query          Actual NCBI query in the NCBI query language. Please
                          refer to the documentation on how to construct one
                          correctly.
    -Q , --query-file     Actual NCBI query in the NCBI query language, stored
                          in a simple text file.
    -r , --rank           Taxonomic rank to perform the identifications on
                          (genus, species, subspecies, serotype) <str> (default:
                          species)
    -s , --sequences      TAB DELIMITED input file containing the the name of
                          the taxon with no special characters, and an
                          underscore '_' instead of spaces, a user defined
                          accession code and the path of the fasta file. The
                          fasta file that the path point to can be either
                          uncompressed or compressed with gzip/bgzip
    -a , --accessions     TAB DELIMITED input file containing the the name of
                          the taxon with no special characters, and an
                          underscore '_' instead of spaces, a user defined valid
                          NCBI nucleotide, assembly or WGS accession code.
    -g  [ ...], --genera  [ ...]
                          List containing the names of specific genera the
                          abundances should be calculated on, separated by a
                          space character <genus1 genus2 genus3 ...>
    -c , --cores          Number of cores for HAYSTACK to use
    -M , --mem            Max memory resources allowed to be used ofr indexing
                          the input for the filtering alignment (default: max
                          available memory 8192.0)
    -u, --unlock          Unlock the working directory after smk is abruptly
                          killed <bool> (default: False)
    -d, --debug           Debug the HAYSTACK workflow <bool> (default: False)
    -smk , --snakemake    Snakemake flags (default: '')

haystac sample
--------------

    -h, --help            show this help message and exit
    -p , --sample-prefix 
                          Sample prefix for all the future analysis. Optional if
                          SRA accession is provided instead <str>
    -o , --sample-output-dir 
                          /path/to/sample <str>
    -f , --fastq          Path to the fastq input file. Can be raw or with
                          adapters removed
    -f1 , --fastq-r1      Path to the mate 1 fastq input file, if reads are PE.
                          Can be raw or with adapters removed
    -f2 , --fastq-r2      Path to the mate 2 fastq input file, if reads are PE.
                          Can be raw or with adapters removed
    -SA , --sra           Fetch raw data files from the SRA using the provided
                          accession code <str>
    -C, --collapse        Collapse paired end reads <bool> (default: False)
    -T , --trim_adapters 
                          Remove adapters from raw fastq files <bool> (default:
                          True)
    -TF , --adaperremoval-flags 
                          Additional flags to provide to Adapterremoval <str>
    -c , --cores          Number of cores for HAYSTACK to use
    -M , --mem            Max memory resources allowed to be used ofr indexing
                          the input for the filtering alignment (default: max
                          available memory 8192.0)
    -u, --unlock          Unlock the working directory after smk is abruptly
                          killed <bool> (default: False)
    -d, --debug           Debug the HAYSTACK workflow <bool> (default: False)
    -smk , --snakemake    Snakemake flags (default: '')
    --dry-run

haystac analyse
---------------

    -h, --help            show this help message and exit
    -m , --mode           Analysis mode for the selected sample
    -D , --database       Path to the database yaml file used for the sample
                          analysis. MANDATORY
    -S , --sample         Path to the sample parameter yaml file. MANDATORY
    -g  [ ...], --genera  [ ...]
                          List containing the names of specific genera the
                          abundances should be calculated on, separated by a
                          space character <genus1 genus2 genus3 ...>
    -o , --analysis-output-dir 
                          Path to results directory.
    -T , --read-probability-threshold 
                          Posterior probability threshold for a read to belong
                          to a certain species. Chose from 0.5, 0.75 and 0.95
                          (default:0.75).
    -c , --cores          Number of cores for HAYSTACK to use
    -M , --mem            Max memory resources allowed to be used ofr indexing
                          the input for the filtering alignment (default: max
                          available memory 8192.0)
    -u, --unlock          Unlock the working directory after smk is abruptly
                          killed <bool> (default: False)
    -d, --debug           Debug the HAYSTACK workflow <bool> (default: False)
    -smk , --snakemake    Snakemake flags (default: '')
    --dry-run


Citations
========= 

bowtie2
mapdamage 
adaperremoval
snakemake 


FAQs 
====


Contributing
============
Evangelos Antonios Dimopoulos,
Evan K. Irving-Pease,
Alberto Carmagnini


Credits
=======

Don't know who to credit


Change Log 
==========


License
=======

MIT


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
