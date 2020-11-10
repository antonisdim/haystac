# Title     : coverage_t_test
# Objective : Function that calculates the pvalue of the coverage. Testing if there was clustering bias during the
#     sequencing/the reads that contribute to the identification/abundance of a species come only
#     from a specific genomic region
# Created by: Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Created on: 10/11/2020

cov_file <- snakemake@input[[1]]
taxon_fasta_idx <- snakemake@input[[2]]
taxon <- snakemake@wildcards[["orgname"]]
outfile <- snakemake@output[[1]]

file_empty <- function(filenames) file.info(filenames)$size == 0

genome_sizes <- function(taxon_fasta_idx) {
  faidx <- read.csv(taxon_fasta_idx, sep='\t')
  names(faidx) <- c('Name','Length','Offset', 'Linebases', 'Linewidth')
  taxon_seq_len <- sum(faidx[,"Length"])
  return(taxon_seq_len)
}

# check if the coverage stats file is empty

if (file_empty(cov_file)) {
  stop(paste("The file with the coverage stats", cov_file, "is empty.", sep = " "))
}

# get the total length of the ref genome
taxon_seqlen <- genome_sizes(taxon_fasta_idx)

# read the cov stats file
cov_stats <- read.csv(cov_file, sep='\t')
names(cov_stats) <- c('observed','expected')

expected_coverage <- cov_stats[1,"expected"]
observed_coverage <- cov_stats[1,"observed"]

contingency_first_row <- c(observed_coverage, expected_coverage)

print(paste("Observed and expected coverage are", contingency_first_row, sep='\t'))

contingency_second_row <- c(taxon_seqlen, taxon_seqlen)

test_input <- matrix(rbind(contingency_first_row, contingency_second_row), nrow=2)

#perform the test
pvalue <- fisher.test(test_input)$p.value

#write to the output file
write.table(data.frame(taxon, pvalue), file=outfile,
            sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

