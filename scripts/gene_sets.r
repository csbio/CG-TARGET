#!/usr/bin/env Rscript

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script quickly prints out all of the gene sets currently
# available to the user.
# It takes one optional argument that narrows the list down to
# only gene sets associated with a particular GI dataset.
library(data.table)

TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/table_printing.r'))

gene_set_config_file = file.path(TARGET_PATH, 'data/gene_sets/gene_set_config.txt')
gene_set_tab = fread(gene_set_config_file, sep = '\t')

dataset = commandArgs(TRUE)[1]
if (!is.na(dataset)) {
    gene_set_tab = gene_set_tab[dataset_name == dataset]
}

# Select only the gene_set name and description columns
gene_set_tab = gene_set_tab[, list(dataset_name, gene_set_name, description)]

# Get the prettified version of the table!
pretty_tab_string = prettify_table(gene_set_tab)

# Write out to the terminal!
write(pretty_tab_string, stdout())



