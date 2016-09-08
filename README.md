#**CG-TARGET**

CG-TARGET is a collection of R scripts to be used for the purpose of predicting gene-set targets from chemical-genetic interaction profiles.

##Features

##Installation

###Requirements

This software is written in R, and thus requires a working R installation. Microsoft Open R (or any other R with an optimized BLAS library, such as OpenBLAS) is recommended, as speed at which some steps finish depends highly on the speed of the matrix multiplications involved.

__**The following libraries are required:**__

hello:

	data.table
	digest
	ggplot2
	grid
	gridExtra
	optparse
	reshape2
	tools
	yaml

__**The following libraries are optional:**__

(they are only required for the `visualize_gene_targets.r` script, which is optional)

ctc
fastcluster

###Downloading CG-TARGET





License
-------
