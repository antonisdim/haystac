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

probabilities: directory where the likelihood matrix and final posterior abundances/probabilities are stored. The final output for abundance calculation has the suffix posterior_abundance.tsv,

mapdamage: directory that includes all the mapdamage profiles for every taxon in our Database

reads: directory including all the Dirichlet reads for each taxon in out database. 