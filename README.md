[![Documentation Status](https://readthedocs.org/projects/haystac/badge/?version=latest)](https://haystac.readthedocs.io/en/latest/?badge=latest)

# HAYSTAC: A Bayesian framework for robust and rapid species identification in high-throughput sequencing data 

## Introduction 

`haystac` is a light-weight, fast, and user-friendly species identification tool. It evaluates the presence of a 
particular species of interest in a metagenomic sample, and provides statistical support for the species assignment. 
The method is designed to estimate the probability that a specific taxon is present in a metagenomic sample given a set 
of sequencing reads and a database of reference genomes. It works equally well with both modern and ancient DNA sequence
data.

## Setup

`haystac` can be run on either macOS or Linux based systems.

The easiest way to install `haystac` and all its dependencies is via the [conda package manager](
https://docs.conda.io/projects/conda/en/latest/index.html).

### Install conda
To install `miniconda3` for macOS:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh
```
or for Linux:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh
```
To install `mamba`"
```bash
conda install mamba -n base -c conda-forge
```
### Install haystac
<!-- Then use `conda` to install `haystac` from the [bioconda](https://bioconda.github.io/) channel: -->
Then use `conda` to install `haystac`:
```
conda install -c bioconda haystac
```
We suggest you install `haystac` into a new environment, as you may encounter install delays if you try to install 
into an existing environment. 

We also advise using `mamba` for the faster installation of `haystac`:
```
mamba install -c bioconda haystac
```

## Quick Start

`haystac` consists of three main modules:
1) `database` for building a database of reference genomes
2) `sample` for pre-processing of samples prior to analysis
3) `analyse` for analysing a sample against a database

### 1. Build a database

To begin using `haystac` we firstly need to construct a database containing all species of interest to our study. In our 
[preprint](https://www.biorxiv.org/content/10.1101/2020.12.16.419085v1), we show that `haystac` makes robust species 
identifications with genus specific databases (for prokaryotes), allowing for very fast hypothesis driven analyses.

In this example, we will build a database containing all species in the *Yersinia* genus, by supplying `haystac` with a 
simple NBCI search query.
```
haystac database \
    --mode build \
    --query '"Yersinia"[Organism] AND "complete genome"[All Fields]' \
    --output yersinia_db
```

To construct an NCBI search query for your area of interest, visit the [NCBI Nucleotide database](
https://www.ncbi.nlm.nih.gov/nucleotide/) and use the search feature to obtain a correctly formatted query string from 
the "Search details" box. This search query can be used directly with `haystac` to automatically download and build 
a reference database based on the accession codes present in the resultset returned by the query.

For more exhaustive analyses, you can build a database containing the 5,681 species present in the [RefSeq 
representative](https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/) database of prokaryotic species by running:
```
haystac database \
    --mode build \
    --refseq-rep prokaryote_rep \
    --output refseq_db
```

Note: Building a database this big is not recommended on a laptop computer. 

### 2. Prepare a sample for analysis

The second step in using `haystac` is to prepare a sample for analysis.

In this example, we will download an aDNA library from [Rasmussen *et al.* (2015)](https://doi.org/10.1016/j.cell.2015.10.009), 
by giving `haystac` the SRA accession code [ERR1018966](https://www.ncbi.nlm.nih.gov/sra/?term=ERR1018966). 
Most published genomics papers include a BioProject code (e.g. [PRJEB10885](
https://www.ncbi.nlm.nih.gov/bioproject/PRJEB10885)), from which you can obtain SRA accessions for each sequencing 
library.
``` 
haystac sample \
    --sra ERR1018966 \
    --output ERR1018966
```

To prepare a sample of your own, you will need either single-end or paired-end short read sequencing data in 
`fastq` format. 

For a paired-end library, you specify the location of the `fastq` files and the name of
the output directory. You may also choose to collapse overlapping mate pairs (e.g. for an aDNA library).
```
haystac sample \
    --fastq-r1 /path/to/sample1_R1.fq.gz \
    --fastq-r2 /path/to/sample1_R2.fq.gz \
    --collapse True \
    --output sample1
```
By default, `haystac` will scan the supplied library, identify adapter sequences, and automatically remove them.

### 3. Analyse a sample against a database

The third step in using `haystac` is to perform an analysis of a sample against a database.

Here, we will use `haystac` to calculate the mean posterior abundance of all species in the *Yersinia* genus found within
the sample `ERR1018966`.
```
haystac analyse \
    --mode abundances \
    --database yersinia_db\
    --sample ERR1018966 \
    --output yersinia_ERR1018966
```

When the analysis is complete, there will be several new sub-folders in the output directory `yersinia_ERR1018966/`. To 
determine if sample `ERR1018966` contains *Yersinia pestis* (i.e. the plague) we can consult the spreadsheet containing
the mean posterior abundance estimates for all species in the *Yersinia* database (i.e., 
`yersinia_ERR1018966/probabilities/ERR1018966/ERR1018966_posterior_abundance.tsv`). From this, we can see that 3,266 
reads were uniquely assigned to *Yersinia pestis*, with an overall abundance of 0.047%, and that the chi-squared test 
indicates that the reads are spread evenly across the genome.

Before we can confidently conclude that `ERR1018966` contains ancient *Yersinia pestis*, we may want to perform a 
damage pattern analysis.
```
haystac analyse \
    --mode abundances \
    --database yersinia_db\
    --sample ERR1018966 \
    --output yersinia_ERR1018966 \
    --mapdamage True
```


## User documentation

`haystac` has many features and potential uses, and we encourage you to use module help menus (e.g. `haystac database --help`) 
to explore these options. The full user documentation is available here: https://haystac.readthedocs.io/en/master/


## Reporting errors

`haystac` is under active development and we encourage you to report any issues you encounter via the [GitHub issue 
tracker](https://github.com/antonisdim/haystac/issues).  

 
## Citation

A preprint describing `haystac` is available on *bioRxiv*:
 
> Dimopoulos, E.A.\*, Carmagnini, A.\*, Velsko, I.M., Warinner, C., Larson, G., Frantz, L.A.F., Irving-Pease, E.K., 2020. 
> HAYSTAC: A Bayesian framework for robust and rapid species identification in high-throughput sequencing data. 
> *bioRxiv* 2020.12.16.419085. https://www.biorxiv.org/content/10.1101/2020.12.16.419085v1

## License 
MIT
