Tutorial
========

Configuring HAYSTAC
-------------------

If a user is running haystac for the first time on a machine, they might want to run the ``haystac config`` module first. The user can use this module to provide an API key for querying NCBI and their preferred path for the cache genomes folder, among other things. All of the options have default values that can be used. If a user later on wishes to change any of these parameters specifically they can either run ``haystac config`` to pass a value to a specific argument.

Here is an example command that allows the configuration of using conda as a package manager for running the other haystac modules.::

    haystac config --use-conda True

Building the database
---------------------

In order to build the database we will be using the ``database`` module of HAYSTAC.
First we need to know what organisms we would like to include in our database. Do we only need the complete genomes of a specific genus or do we want more genera? 

Constructing the Query
----------------------

After deciding what taxa we would like to include in our database, we need to construct an NCBI query that will return all the accessions that belong to the taxa that we are interested in. One of the best ways to construct such a query is to go on the website of  NCBI's Nucleotide database (https://www.ncbi.nlm.nih.gov/nucleotide/), type in our query and get the correctly formatted query string from the "Search details" box. We can then use that string for the construction of our database. 

For example if we would like to build a database of all the complete genomes of the species in the Yersinia genus we can use the following command:::

    haystac database --mode build \
        --query '"Yersinia"[Organism] AND "complete genome"[All Fields]' \
        --output yersinia_example

For each species (or any other user defined taxonomic rank), the longest sequence per taxon will be used to populate our database. 

Representative RefSeq species
-----------------------------

When constructing a database there is always the option to include the species of the prokaryotic representative RefSeq database as well. All you need to do is include the corresponding flag in your command.::

    haystac database --mode build \
        --query '"Yersinia"[Organism] AND "complete genome"[All Fields]' \
        --output yersinia_example \
        --refseq-rep prokaryote_rep

Important note on RefSeq databases
----------------------------------

``haystac database`` currently can build databases from three RefSeq tables, the prokaryotic representative RefSeq table, the eukaryotes RefSeq table and the viruses RefSeq table. When the prokaryotic representative the database is built, only species of microorganisms are included (strains are excluded), whereas in the eukaryotes and viruses databases subspecies and strains are included respectively. To build any of the above databases, specify the desired RefSeq table to be used by the ``--refseq-rep`` flag (``prokaryote_rep`` for the prokaryotic representative, ``eukaryotes`` for the eukaryotes and ``viruses`` for the viruses table).

Providing custom accessions
---------------------------

It is also possible to provide your own accessions for a selected species/taxon. For that you will need to prepare a tab delimited file with two columns. The first column is the name of the taxon, that cannot contain any special characters, other than an underscore ('_'), and the second column is a valid NCBI accession code. 

Here is an example of the contents of such a file:::

    Yersinia_ruckeri    NZ_CP025800.1

The we can simply run the following command:::

    haystac database --mode build \
        --query '"Yersinia"[Organism] AND "complete genome"[All Fields]' \
        --output yersinia_example \
        --accessions-file acc_example.txt

Providing custom sequences
--------------------------

It is also possible to provide your own sequences for a taxon. To do that you will need a a tab delimited file containing the the name of the taxon with no special characters, and an underscore ('_') instead of spaces, a user defined accession code and the path of the fasta file. The fasta file that the path point to can be either uncompressed or compressed with gzip/bgzip.

Here is an example of such a file:::

    Yersinia_ruckeri    user_seq_1  ~/example_sequence/user_seq_1.fasta

The we can simply run the following command:::

    haystac database --mode build \
        --query '"Yersinia"[Organism] AND "complete genome"[All Fields]' \
        --db-output yersinia_example \
        --sequences-file seq_example.txt

Combinations
------------

All of the previous options can be combined into one command. It is important to note that only one sequence file per taxon is allowed in our database, and the priority goes user defined accessions or sequences, representative RefSeq and then user specified query. It should be noted that for custom NCBI queries, plasmids can be fetched if they are part of a genome assembly. The only exception to the one sequence file per taxon rule is plasmid sequences of the RefSeq representative that are not part of complete genome assemblies and for that reason they are downloaded separately. 

Index building 
--------------

For the first part of the analysis an index out of all the genomes that are included in our database needs to be build. This is a process that can take big amounts of memory depending on the number and the complexity of sequences that our database includes. For that reason the user can specify the desired amount of memory resources available to ``haystac`` and the program will try to build the required index. This can be specified through the ``--mem`` flag, that can be appended to the any of the commands shown above. Memory resources need to be specified in MB. If the memory resources provided are less than the size of the files that need to be indexed an error will be raised. We also advise caution when changing the bowtie2 file size scaling factor.

Database building modes
-----------------------

For the complete construction of a database, sequences need to be downloaded and subsequently indexed.
By specifying ``--mode build`` to ``haystac database``, the program downloads and indexes all the sequences that have been requested by the user in one step.
If a user would like to only download sequence data and index them later it is possible to do so, by specifying ``haystac database --mode fetch``, to download the sequences first and then execute ``haystac database --mode index`` in order to perform the indexing.
If mode ``fetch`` is run first then mode ``index`` should be run subsequently, and not mode ``build``, otherwise an error will be raised.

Here is an example of building a database in two steps instead of one:::

    haystac database --mode fetch \
        --query '"Yersinia"[Organism] AND "complete genome"[All Fields]' \
        --output yersinia_example
    haystac database --mode index \
        --output yersinia_example

Building a mitochondrial DNA database
-------------------------------------

When a user is providing a query about eukaryotes it is also possible to build a database with only mitochondrial genomes (by default whole genome assemblies will be fetched for a given query). In order to do that a user can specify the ``--mtDNA`` flag when running ``haystac database``. We strongly advise against having a mixed database of full eukaryotic genome assemblies for certain taxa and only mtDNA sequences for other taxa, as this will bias the identifications towards the taxa with full genome assemblies.

Preparing a sample for analysis
-------------------------------

After our database is built we need to prepare our samples for analysis. For that purpose, we are using the ``sample`` module of haystac. The input files can be SE, PE or collapsed reads. If the reads are collapsed they are going to be treated as SE reads.

It is possible to trim sequencing adapters and collapse PE reads by specifying the relative flags. Samples (specific sequencing runs) can be also downloaded from the SRA if an sra run accession is provided. 

If you have SE or already collapsed reads you only need to specify a file path for the ``--fastq`` flag.
If your input is PE reads then you will need to specify file paths for both the ``--fastq-r1`` and ``--fastq-r2``.
If you want to download files from the SRA all you need to do is provide an SRA accession for the ``--sra`` flag.

Here is an example of downloading reads from the SRA, trimming sequencing adapters and collapsing reads.::

    haystac sample --sra ERR1018966 \
        --output sample_example

Sample analysis
---------------

In order to analyse any sample we will need to use the ``analyse`` module of haystac.

Filtering Alignment
-------------------

The first step for the sample analysis is to filter in all the reads that align to any of the genomes in our database. For that we will need to use the ``haystac analyse --mode filter``.

Here is an example command:::

    haystac analyse --mode filter \
        --database yersinia_example \
        --sample sample_example \
        --output analysis_output

Database Alignments
-------------------

After we have filtered our libraries we can align the filtered reads against all the genomes that are included in our database. This can be done by using mode ``align`` of ``haystac analyse``.

For example:::

    haystac analyse --mode align \
        --database yersinia_example \
        --sample sample_example \
        --output analysis_output

Unless the user has a deep understanding of their dataset we advise to be cautious when changing the base mismatch probability that is used later on in the method's probabilistic model.

Likelihood calculation
----------------------

After all the individual alignments have been competed, the number of transitions and transversions will be counted for every read that has aligned against any of the reference genomes in our database. Then the likelihoods and posterior probabilities for each read being sampled from a given reference genome will be calculated. For this step we can use the ``likelihoods`` mode of ``haystac analyse``.::

    haystac analyse --mode likelihoods \
        --database yersinia_example \
        --sample sample_example \
        --output analysis_output

Important Note on the Dirichlet Assignment process during Likelihood calculation
--------------------------------------------------------------------------------

It is important to be aware of the individual read posterior probability threshold, for a read to be assigned to a taxon. As a default HAYSTAC uses the conservative 0.75 probability threshold for the Dirichlet assignment. The higher value you pick the more conservative the assignments become. It is useful to sometimes pick a value depending on what taxa are being identified. If there is a need to distinguish between closely related taxa then a more conservative threshold would increase the specificity of the analysis therefore being more appropriate, whereas when you're trying to generally characterise a metagenome a less conservative value could increase the sensitivity of the analysis be more helpful.

Single organism sample or metagenome ? 
--------------------------------------

Depending on whether we would like to identify the species a sample is belongs to, or perform a metagenomic analysis, we can use the ``probabilities`` or ``abundances`` mode of ``haystac analyse`` respectively.

Assignment Probability Calculation
----------------------------------

In order to calculate posterior assignment probabilities we can run the following command:::

    haystac analyse --mode probabilities \
        --database yersinia_example \
        --sample sample_example \
        --output analysis_output

Mean Posterior Abundances
-------------------------

In order to calculate mean posterior abundances we can run the following command:::

    haystac analyse --mode abundances \
        --database yersinia_example \
        --sample sample_example \
        --output analysis_output

Along with the abundance calculation, we also perform a chi2 test to assess if the reads that have been assigned to a taxon are clustering around specific genomic areas or if they represent a random sample of the organism's genome. The results of this test should be trusted for low depth sequencing data (equal or less than 1X). The null hypothesis is that there is no read clustering.

Reads
-----

After the mean posterior abundances have been calculated for a sample, all the reads that have been assigned to a taxon through the Dirichlet process can be outputted in separate bam files ready for further downstream analyses (like assembling or variant calling for instance) via the ``reads`` module. Reads that have been assigned to the Grey and Dark Matter are outputted in fastq files as they have not been uniquely assigned to a taxon.

Here is an example command:::

    haystac analyse --mode reads \
        --database yersinia_example \
        --sample sample_example \
        --output analysis_output

Mapdamage analysis
------------------

If our samples are ancient we can use mapDamage to estimate the level of deamination in the reads that have aligned to any taxon in our database. For that we can use the ``mapdamage`` module of haystac. The mapDamage analysis will be performed on the subset of reads that have been uniquely assigned to a taxon through the dirichlet process. This module can be either run independently or after the ``reads`` module.

Here is an example command:::

    haystac analyse --mode mapdamage \
        --database yersinia_example \
        --sample sample_example \
        --output analysis_output

Important note on sample analysis
---------------------------------

The first 3 steps (modes: ``filter``, ``align``, ``likelihoods``) can be executed automatically when you call the ``probabilities`` or ``abundances`` mode of haystac.


