#!/usr/bin/env Rscript
# This script takes in a gene-target prediction table and
# exports a clustergram of the gene-target prediction
# scores. Row and column labels are specified using
# arguments for row/column info tables and the columns
# that I want to extract from these tables.

library(getopt)
library(reshape2)
library(data.table)
library(fastcluster)
library(ctc)

# Source in the similarity metric(s) I will use 
# NOTE: I am not computing p values
TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/cosine.r'))
source(file.path(TARGET_PATH, 'lib/dist_for_clustering.r'))
source(file.path(TARGET_PATH, 'lib/filenames.r'))

spec = rbind(c('dataset_file', 'f', 1, 'character', 'A 4-column, gzipped table, with one column for the row identifiers (query_key), two columns for the column identifiers (screen_name, expt_id), and one column with gene-target prediction scores (score).'),
             c('clustergram_name', 'n', 1, 'character', 'The name of the cdt/atr/gtr files that are exported.'),
             c('GI_info_table', 'G', 1, 'character', 'A table, with one column as a unique query identifier (query_key) and the remaining columns providing information linked to that query identifier (systematic gene name, interpretable gene name, etc.'),
             c('sample_table', 'S', 1, 'character', 'The sample table used to process the chemical genomic data from the beginning. Must have "screen_name" and "expt_id" columns.'),
             c('GI_columns', 'g', 1, 'character', 'Columns from the GI_info_table to include in the clustergram row labels. Comma-delimited, no space between.'),
             c('condition_columns', 'c', 1, 'character', 'Columns from the sample_table to include in the clustergram column labels. Comma-delimited, no space between.'),
             c('verbosity', 'v', 2, 'character', 'Columns from the sample_table to include in the clustergram column labels. Comma-delimited, no space between.')
             )

opt = getopt(spec)

# Load in the dataset
dataset_file = normalizePath(opt$dataset_file)
tab = fread(sprintf('gzip -dc %s', dataset_file), header = TRUE, colClasses = c('character', 'character', 'character', 'numeric'))

# Define output folder
output_folder = sub('.txt.gz', '', dataset_file)
dir.create(output_folder)

# Read in GI info and sample tables
gi_info_tab = fread(opt$GI_info_table, header = TRUE, colClasses = 'character')
setkey(gi_info_tab, query_key)
print(gi_info_tab)

sample_tab = fread(opt$sample_table, header = TRUE, colClasses = 'character')
sample_tab[, condition_key := sprintf('%s_%s', screen_name, expt_id)]
setkey(sample_tab, condition_key)
print(sample_tab)

# Reshape the table into a matrix, doing appropriate things to the row/col
# names to keep things lined up.
tab[, condition_key := sprintf('%s_%s', screen_name, expt_id)]
mat = acast(data = tab, formula = query_key ~ condition_key, value.var = 'score')

# Get the final rownames of the matrix!
GI_columns = strsplit(opt$GI_columns, ',')[[1]]
print(GI_columns)
row_labels = apply(gi_info_tab[dimnames(mat)[[1]]][, GI_columns, with = FALSE], 1, paste, collapse = '_')

# Get the final colnames of the matrix!
condition_columns = strsplit(opt$condition_columns, ',')[[1]]
print(condition_columns)
col_labels = apply(sample_tab[dimnames(mat)[[2]]][, condition_columns, with = FALSE], 1, paste, collapse = '_')

# Change the dimnames of the matrix!
dimnames(mat) = list(row_labels, col_labels)

# Compute distance matrices and cluster!
dists = list(rows = as.dist(cos_dist_for_clust(mat, na.rm = TRUE)),
             cols = as.dist(cos_dist_for_clust(t(mat), na.rm = TRUE))
             )

hclusts <- lapply(dists, hclust, method = 'average')
names(hclusts) <- names(dists)

# Get the output filename
near_full_path = file.path(output_folder, opt$clustergram_name)

r2gtr(hr = hclusts$rows, file = sprintf('%s.gtr', near_full_path), distance = 'custom')
r2atr(hc = hclusts$cols, file = sprintf('%s.atr', near_full_path), distance = 'custom')
r2cdt(hr = hclusts$rows, hc = hclusts$cols, data = mat, file = sprintf('%s.cdt', near_full_path))


