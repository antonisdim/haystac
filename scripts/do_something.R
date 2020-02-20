#!/usr/bin/env Rscript

suppressMessages(library(rjson, quietly=T, warn.conflicts=F))

# Do something
#
do_something <- function(input_file, output_file) {
	blah = c()

    cat(toJSON(blah, indent = 2), file=output_file)
}


if (sys.nframe() == 0) {
	# redirect all output to the log
	log <- file(snakemake@log[[1]], open="wt")
	sink(log)

	do_something(
		input_file=snakemake@input[[1]],
		output_file=snakemake@output[[1]]
	)
}
