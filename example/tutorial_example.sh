#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# installing haystac
python -m pip install git+https://github.com/antonisdim/haystac.git

# configure haystac by providing a pth for the cache folder
haystac config --cache ./haystac_cache_example/

# create an example database for the genus yersinia
haystac database --mode build \
    --query '("Yersinia"[Organism] OR "Yersinia"[Organism]) AND "complete genome"[All Fields]' \
    --output yersinia_example

# download a sample from the SRA and prepare it for analysis
haystac sample --sra ERR1018966 --output sample_example

# perform a filtering alignment against our Yersinia database
haystac analyse --mode filter \
    --database yersinia_example \
    --sample sample_example \
    --output analysis_output

# perform the individual metagenomic alignments
haystac analyse \
    --mode align \
    --database yersinia_example \
    --sample sample_example \
    --output analysis_output

# perform the likelihood calculation
haystac analyse --mode likelihoods \
    --database yersinia_example \
    --sample sample_example \
    --output analysis_output

# calculate the posterior abundances of the Yersinia species in our sample
haystac analyse --mode abundances \
    --database yersinia_example \
    --sample sample_example \
    --output analysis_output

# output all the assigned reads in individual bam files
haystac analyse --mode reads \
    --database yersinia_example \
    --sample sample_example \
    --output analysis_output

# perform a mapdamage analysis for each taxon in out database
haystac analyse --mode mapdamage \
    --database yersinia_example \
    --sample sample_example \
    --output analysis_output