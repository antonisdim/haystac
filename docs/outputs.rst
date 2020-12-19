Outputs
=======

This is a list of all the outputs produced after a successful run of each of ``haystac``'s modules. All these outputs will be found under the ``--output`` path specified by the user when running any of ``haystac``'s modules. We advise to create separate output directories for each one of the modules (``database``, ``sample`` and ``analyse``).

Expected outputs for ``haystac database``
-----------------------------------------

**db_taxa_accessions.tsv**: includes all the taxon/accession pairs that are included in a database.

**idx_database.done**: file indicating that all the individual bowtie2 indices for each taxon have been prepared

**entrez**: directory containing all the results from querying the NCBI, including the nucleotide and taxonomy databases

**bowtie**: directory containing the bowtie2 index files for all out database, that will be used for the filtering alignment

**database_inputs**: directory containing the representative RefSeq table that is downloaded from NCBI.

Expected outputs for ``haystac sample``
---------------------------------------

**fastq_inputs**: folder containing the outputs of the sample module.

**fastq_inputs/meta**: directory that includes the read count file.

**fastq_inputs/(SE | PE | COLLAPSED)**: directory containing the trimmed reads produced by AdapteRremoval

Expected outputs for ``haystac analyse``
----------------------------------------

**bam**: directory containing the bam file of the filtering alignment

**fastq**: directory containing the filtered reads in fastq format along with their average read length

**alignments**: directory where all the individual alignment bam files for each taxon in a database of the filtered reads are outputted

**ts_tv_counts**: directory where all the transition and transversion counts are stored per taxon

**probabilities**: directory where the likelihood matrix and final posterior abundances/probabilities are stored. The final output for abundance calculation has the suffix posterior_abundance.tsv,

**mapdamage**: directory that includes all the mapdamage profiles for every taxon in our Database

**reads**: directory including all the Dirichlet reads for each taxon in out database.