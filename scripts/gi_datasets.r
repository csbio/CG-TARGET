#!/usr/bin/env Rscript

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script quickly prints out all of the genetic interaction
# datasets currently available (or downloadable) for the user.
library(data.table)

TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/table_printing.r'))

dataset_config_file = file.path(TARGET_PATH, 'data/GI/dataset_config.txt')
dataset_tab = fread(dataset_config_file, sep = '\t')

# Select only the dataset name and description columns
dataset_tab = dataset_tab[, list(dataset_name, description)]

# Get the prettified version of the table!
pretty_tab_string = prettify_table(dataset_tab)

# Write out to the terminal!
write(pretty_tab_string, stdout())
