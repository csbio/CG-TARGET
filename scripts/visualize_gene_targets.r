#!/usr/bin/env Rscript

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script takes in a gene-target prediction table and
# exports a clustergram of the gene-target prediction
# scores. Row and column labels are specified using
# arguments for row/column info tables and the columns
# that I want to extract from these tables.

library(optparse)
library(yaml)
library(reshape2)
library(data.table)
library(fastcluster)
library(ctc)

TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/cosine.r'))
source(file.path(TARGET_PATH, 'lib/dist_for_clustering.r'))
source(file.path(TARGET_PATH, 'lib/filenames.r'))
source(file.path(TARGET_PATH, 'lib/pos_args.r'))
source(file.path(TARGET_PATH, 'lib/datasets.r'))

positional_arguments_list = c(CONFIG_FILE = 'yaml-formatted gene-set prediction configuration file.')

option_list = list(make_option('--rand', action = 'store_true', default = FALSE, help = 'Use this flag if you want to perform visualization for the predicted gene-level targets of the resampled profiles')
                                      )

final_usage_positional = get_usage_and_positional(positional_arguments_list)
parser = OptionParser(usage = final_usage_positional, option_list = option_list)
arguments = parse_args(parser, positional_arguments = length(positional_arguments_list))
print(arguments)
arg = arguments$args
opt = arguments$options


# Read in the configuration parameters
config_f = file(arg[1], 'rt')
config_params = yaml.load_file(config_f)
close(config_f)

output_folder = config_params$Required_arguments$output_folder

gene_target_dir = get_gene_target_folder(output_folder)
if (opt$rand) {
    gene_target_filename = get_gene_target_prediction_resampled_filename(
                                gene_target_dir,
                                config_params$Required_arguments$`per-array_resampling_scheme`,
                                config_params$Required_arguments$`num_per-array_resampled_profiles`,
                                config_params$Required_arguments$`per-array_resampling_seed`,
                                config_params$Required_arguments$cg_gi_similarity_measure
                                )
} else {
    gene_target_filename = get_gene_target_prediction_filename(
                                                       gene_target_dir,
                                                       config_params$Required_arguments$cg_gi_similarity_measure
                                                       )
}


# Load in the dataset
dataset_file = normalizePath(gene_target_filename)
tab = fread(sprintf('gzip -dc %s', dataset_file), header = TRUE, colClasses = c('character', 'character', 'character', 'numeric'))

# Read in GI info
gi_info = get_gi_info(config_params$Required_arguments$gi_dataset_name, TARGET_PATH)
#gi_tab = fread(sprintf('gzip -dc %s', gi_filenames$gi_tab))
query_info_tab = fread(gi_info$gi_query_tab, header = TRUE, colClasses = 'character')

setkey(query_info_tab, query_key)
print(query_info_tab)

# Read in the CG info
if (opt$rand) {
    # If this for the resampled profiles, which do not have an accompanying
    # sample table, then just make one on the fly!
    sample_tab = unique(tab[, list(screen_name, expt_id)])
} else {
    # Otherwise read in the real sample table
    sample_tab = fread(config_params$Required_arguments$cg_col_info_table, header = TRUE, colClasses = 'character')
}

sample_tab[, condition_key := sprintf('%s_%s', screen_name, expt_id)]
setkey(sample_tab, condition_key)
print(sample_tab)

# Reshape the table into a matrix, doing appropriate things to the row/col
# names to keep things lined up.
tab[, condition_key := sprintf('%s_%s', screen_name, expt_id)]
mat = acast(data = tab, formula = query_key ~ condition_key, value.var = 'score')

# Get the final rownames of the matrix!
GI_columns = strsplit(config_params$Options$gene_target_prediction_visualization$gi_col_table_vis_columns, ',')[[1]]
print(GI_columns)
row_labels = apply(query_info_tab[dimnames(mat)[[1]]][, GI_columns, with = FALSE], 1, paste, collapse = '_')

# Get the final colnames of the matrix!
#message(opt$rand)
if (opt$rand) {
    condition_columns = c('screen_name', 'expt_id')
} else {
    condition_columns = strsplit(config_params$Options$gene_target_prediction_visualization$cg_col_table_vis_columns, ',')[[1]]
}
print(condition_columns)
col_labels = apply(sample_tab[dimnames(mat)[[2]]][, condition_columns, with = FALSE], 1, paste, collapse = '_')

#message(str(mat))
#message(str(row_labels))
#message(str(col_labels))

# Change the dimnames of the matrix!
dimnames(mat) = list(row_labels, col_labels)

# Compute distance matrices and cluster!
dists = list(rows = as.dist(cos_dist_for_clust(mat, na.rm = TRUE)),
             cols = as.dist(cos_dist_for_clust(t(mat), na.rm = TRUE))
             )

hclusts <- lapply(dists, hclust, method = 'average')
names(hclusts) <- names(dists)

# Get the output filename
cdt_folder = sub('.txt.gz', '', gene_target_filename)
dir.create(cdt_folder)
near_full_path = file.path(cdt_folder, basename(cdt_folder))

r2gtr(hr = hclusts$rows, file = sprintf('%s.gtr', near_full_path), distance = 'custom')
r2atr(hc = hclusts$cols, file = sprintf('%s.atr', near_full_path), distance = 'custom')
r2cdt(hr = hclusts$rows, hc = hclusts$cols, data = mat, file = sprintf('%s.cdt', near_full_path))