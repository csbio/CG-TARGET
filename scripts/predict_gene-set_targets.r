#!/usr/bin/env Rscript
## CHANGE THIS HEADING DOCUMENTATION!!!!!!!!
#
#
## This script will take in the results of chemical genomics target prediction
## and return "z scores" describing the relationship between a drug's average
## target prediction score for genes that match a GO term, and the drug's overall
## average target prediction score.
#
## This script exports all zscore calculations
## Analysis is done using other scripts
#


library(optparse)
library(yaml)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)

# Source in libraries for gene set target prediction
TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/go_drivers.r'))
source(file.path(TARGET_PATH, 'lib/go_matrices.r'))
source(file.path(TARGET_PATH, 'lib/file_naming.r'))
source(file.path(TARGET_PATH, 'lib/top_table.r'))
source(file.path(TARGET_PATH, 'lib/table_to_globular.r'))
source(file.path(TARGET_PATH, 'lib/zscore-pval_by-drug-and-goterm_fast.r'))
source(file.path(TARGET_PATH, 'lib/pos_args.r'))
source(file.path(TARGET_PATH, 'lib/filenames.r'))

positional_arguments_list = c(CONFIG_FILE = 'yaml-formatted gene-set prediction configuration file.')

option_list = list(make_option('--test_run', action = 'store_true', default = FALSE, help = 'Performs predictions for a small subset of the data, to ensure that everything is working properly before the full set of predictions are made (can take days for large screens). The results from this run will not be accurate!!!')
                   )

#positional_arguments_list = c(gene_predictions = 'Gzipped file containing the gene-level target predictions for each condition.',
#                              rand_gene_predictions = 'Gzipped file containing the gene-level target predictions for each resampled condition.',
#                              gene_sets = 'Two-column, tab-delimited file (WITH HEADER) that maps gene identifiers (must have same column name as corresponding column in the query info table) to gene sets.',
#                              query_info_table = 'Table with a "query_key" column that maps to the query identifiers in the genetic interaction data and another column that maps to the gene identifiers in the gene_sets file.',
#                              driver_cutoff = 'Similarity value (similarity to the condition profile) above which a gene is reported as a "driver" of particular gene-set prediction. Default value cannot be specified because unbounded similarity measures are accepted as input.',
#                              output_table = 'File into which the results are written. (it is gzipped, so the extension should be *.txt.gz)'
#                              )

#option_list = list(
#    make_option('--sample_table', help = 'File with a table that maps conditions to control information'),
#    make_option('--neg_control_column', help = 'Name of True/False column in the sample_table that indicates which samples are to be used as negative controls in gene-set prediction. If this option is not specified, only the resampled conditions will be used for gene-set prediction.'),
#    make_option('--condition_name_column', help = 'Name of a column in the sample table with human-readable condition names'),
#    make_option(c('-n', '--per_condition_randomizations'), type = 'integer', help = 'The number of "per-condition" randomizations to perform. When not specified, the default is to perform the same number of randomizations that are present in the rand_gene_predictions file.'),
#    make_option('--gene_set_names', help = 'Two-column, tab-delimited file that maps each gene set identifier to an interpretable name for gene-set prediciton output. The values in the first column must match the gene sets in the gene_sets file above. The values in the second column will be included in the final results table (along with those in the first column).'),
#    make_option('--gene_name_column', help = 'Column in "query_info_table" that contains interpretable gene names. In S. cerevisiae, for example this would be the column that contains "TUB3" in the row for ORF YML124C'),
#    make_option('--test_run', action = 'store_true', default = FALSE, help = 'Performs predictions for a small subset of the data, to ensure that everything is working properly before the full set of predictions are made (can take days for large screens). The results from this run will not be accurate!!!'),
#    make_option('--num_cores', type = 'integer', default = 1, help = 'The number of cores used to run the gene set predictions!'),
#    make_option(c('-v', '--verbosity'), type = 'integer', default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1, is default.')
#    )

final_usage_positional = get_usage_and_positional(positional_arguments_list)

parser = OptionParser(usage = final_usage_positional, option_list = option_list)
arguments = parse_args(parser, positional_arguments = length(positional_arguments_list))
opt = arguments$options
arg = arguments$args

print(opt)
print(arg)

config_f = file(arg[1], 'rt')
config_params = yaml.load_file(config_f)
close(config_f)

################# Load in the required files, deal with arguments
output_folder = config_params$Required_arguments$output_folder
gene_pred_folder = get_gene_target_folder(output_folder)
gene_pred_file = get_gene_target_prediction_filename(gene_pred_folder,
                                                     config_params$Required_arguments$cg_gi_similarity_measure
                                                     )
gene_prediction_tab = fread(sprintf('gzip -dc %s', gene_pred_file), header = TRUE, colClasses = c('character', 'character', 'character', 'numeric'))

rand_gene_pred_file = get_gene_target_prediction_resampled_filename(
                          gene_pred_folder,
                          config_params$Required_arguments$`per-array_resampling_scheme`,
                          config_params$Required_arguments$`num_per-array_resampled_profiles`,
                          config_params$Required_arguments$`per-array_resampling_seed`,
                          config_params$Required_arguments$cg_gi_similarity_measure
                          )
rand_gene_prediction_tab = fread(sprintf('gzip -dc %s', rand_gene_pred_file), header = TRUE, colClasses = c('character', 'character', 'character', 'numeric'))

sample_table_filename = config_params$Required_arguments$cg_col_info_table
sample_table = fread(sample_table_filename, header = TRUE, colClasses = 'character')

gene_set_file = config_params$Required_arguments$gene_set_table
gene_set_tab = fread(gene_set_file, header = TRUE, colClasses = 'character')

query_info_file = config_params$Required_arguments$gi_query_info_table
query_info_tab = fread(query_info_file, header = TRUE, colClasses = 'character')

driver_cutoff = as.numeric(config_params$Required_arguments$driver_cutoff)
if (is.na(driver_cutoff)) { stop('driver_cutoff argument must be numeric.') }

output_table_folder = get_gene_set_target_prediction_folder(output_folder)
output_table = get_gene_set_target_prediction_filename(
                   output_table_folder,
                   config_params$Required_arguments$`per-array_resampling_scheme`,
                   config_params$Required_arguments$`num_per-array_resampled_profiles`,
                   config_params$Required_arguments$`per-array_resampling_seed`,
                   config_params$Required_arguments$cg_gi_similarity_measure,
                   config_params$Required_arguments$`per-condition_randomization_seed`,
                   config_params$Required_arguments$`num_per-condition_randomizations`,
                   config_params$Required_arguments$gene_set_name,
                   opt$test_run
                   )

dir.create(output_table_folder, recursive = TRUE)

true_false_control_map = c('TRUE' = 'expt_control', 'FALSE' = 'treatment', 'True' = 'expt_control', 'False' = 'treatment')

################### Load in optional files and options
# Gene set ID to interpretable name table
if (!is.null(config_params$Options$gene_set_target_prediction$gene_set_name_table)) {
    gene_set_name_file = config_params$Options$gene_set_target_prediction$gene_set_name_table
    gene_set_name_tab = fread(gene_set_name_file, header = FALSE)
    print(gene_set_name_tab)
    setnames(gene_set_name_tab, c('gene_set', 'gene_set_name'))
} else {
    gene_set_name_tab = NULL
}

if (!is.null(config_params$Options$gene_set_target_prediction$query_genename_column)) {
    gene_name_column = config_params$Options$gene_set_target_prediction$query_genename_column
} else {
    gene_name_column = NULL
}

################ Read in the sample table and negative control/condition name columns, if they exist
# Handle the negative control column
neg_control_col = config_params$Options$gene_set_target_prediction$negative_control_column
if (!is.null(neg_control_col)) {
    if (! neg_control_col %in% names(sample_table)) {
        stop(sprintf('column %s not in cg_col_info_table!', neg_control_col))
    }
    control_map = true_false_control_map[sample_table[[neg_control_col]]]
    names(control_map) = sample_table[, sprintf('%s_%s', screen_name, expt_id)]
} else {
    control_map = rep('treatment', dim(sample_table)[1])
    names(control_map) = sample_table[, sprintf('%s_%s', screen_name, expt_id)]
}

# Handle the condition name column
cond_name_col = config_params$Options$gene_set_target_prediction$condition_name_column
if (!is.null(cond_name_col)) {
    ### May need to revisit if the columns are split or not...
    if (! cond_name_col %in% names(sample_table)) {
        stop(sprintf('column %s not in cg_col_info_table!', cond_name_col))
    }
    condition_name_tab = sample_table[, list(condition = sprintf('%s_%s', screen_name, expt_id), condition_name = sample_table[[cond_name_col]])]
    
    # The condition_tab needs to have rand-by-strain conditions added so the join
    # doesn't cause use to lose all or the rand-by-strain conditions!
    rand_condition_name_tab = unique(rand_gene_prediction_tab[, list(screen_name, expt_id)])[, list(condition = sprintf('%s_%s', screen_name, expt_id), condition_name = sprintf('%s_%s', screen_name, expt_id))]
    condition_name_tab = rbind(condition_name_tab, rand_condition_name_tab)
    print(condition_name_tab)
} else {
    condition_name_tab = NULL
}

print(sample_table)
print(control_map)

################### Smashing some data around
all_prediction_tab = rbind(gene_prediction_tab, rand_gene_prediction_tab)
all_prediction_mat = acast(all_prediction_tab, screen_name + expt_id ~ query_key, value.var = 'score')

################### Line up the control vectors to the matrix with all predictions
rand_control_map = rep('rand-by-strain', dim(unique(rand_gene_prediction_tab[, list(screen_name, expt_id)]))[1])
names(rand_control_map) = unique(rand_gene_prediction_tab[, sprintf('%s_%s', screen_name, expt_id)])

# Check for real and random conditions overlapping (should NEVER happen)
print(control_map[1:10])
print(rand_control_map[1:10])
overlap_conds = intersect(names(control_map), names(rand_control_map))
if (length(overlap_conds) > 0) {
    stop(sprintf('The following real and random conditions have the same name:\n%s', paste(overlap_conds, collapse = '\n')))
}

# Combine the control maps and get one master control type vector!!!
all_control_map = c(control_map, rand_control_map)
all_control_vec = all_control_map[rownames(all_prediction_mat)]
print(all_control_vec)
print(sprintf('Number of %s conditions: %s', names(table(all_control_vec)), table(all_control_vec)))

# If any sample types are still NA, filter them out of the dataset
inds_to_remove = is.na(all_control_vec)
conds_to_remove = rownames(all_prediction_mat)[inds_to_remove]

all_control_vec = all_control_vec[!is.na(all_control_vec)]
all_prediction_mat = all_prediction_mat[!is.na(all_control_vec), ]

# Get a matrix with one column per unique query gene (there may be multiple alleles for the same gene),
# and for each row (condition) in that column, we take the maximum prediction score against all
# alleles of the same query gene. Also, get a matrix that indicates which allele the selected
# max score came from.

# First, get the vector of query genes that define where the multiple-allele cases are.
query_column = names(gene_set_tab)[1]
query_map = query_info_tab[[query_column]]
names(query_map) = query_info_tab[['query_key']]
query_gene_vec = query_map[colnames(all_prediction_mat)]

print(str(all_prediction_mat))
print(str(query_gene_vec))
print(all_prediction_mat[1:5, (dim(all_prediction_mat)[2]-9):(dim(all_prediction_mat)[2])])
print(query_gene_vec[(length(query_gene_vec) - 9):length(query_gene_vec)])

# Get the matrix with the best predictions among multiple alleles
best_query_mat_list = best_score_per_col_group(all_prediction_mat, query_gene_vec)

best_prediction_mat = best_query_mat_list[['best_scores']]
best_query_mat = best_query_mat_list[['best_queries']]

# Get a map from query key column to an interpretable name, if the
# gene_name_column exists. Then create a matrix that matches the
# best_query_mat matrix and contains interpretable query gene
# names. If that column was not given, create a matrix of NA
# values instead.
if (!is.null(gene_name_column)) {
    if (! gene_name_column %in% names(query_info_tab)) {
        stop(sprintf('column %s not in query_info_tab!', gene_name_column))
    }
    gene_name_map = query_info_tab[[gene_name_column]]
    names(gene_name_map) = query_info_tab[['query_key']]
    best_query_name_mat = matrix(gene_name_map[best_query_mat], nrow = nrow(best_query_mat), ncol = ncol(best_query_mat), dimnames = dimnames(best_query_mat))
} else {
    best_query_name_mat = matrix(NA_character_, nrow = nrow(best_query_mat), ncol = ncol(best_query_mat), dimnames = dimnames(best_query_mat))
}

###### Determine which sample types exist and will be used to make predictions
sample_types_to_predict_for <- unique(all_control_vec)

#########################
#########################

#########
### What are the control types? (There may be experimental controls, or maybe not!)
controls <- intersect(sample_types_to_predict_for, c('expt_control', 'rand-by-strain'))

# Set up parallelization
ncores = as.numeric(config_params$Options$gene_set_target_prediction$num_cores)
if (length(ncores) == 0) {
    ncores = 1L
    registerDoSEQ()
} else if (ncores > 1) {
    registerDoParallel(cores = ncores)
} else {
    stop('option num_cores could not be converted into an integer!')
}

# Reshape the gene set table into a matrix
gene_set_mat_formula = as.formula(sprintf('%s ~ %s', names(gene_set_tab)[1], names(gene_set_tab)[2]))
gene_set_tab[, annotation := 1]
gene_set_matrix = acast(data = gene_set_tab, formula = gene_set_mat_formula, value.var = 'annotation', fill = 0)

## Get a matrix of all possible propagated orf to GO term mappings
#raw_go_matrix <- get_gene_go_bp_is_a_propagated_matrix(raw_orfs)
#
## Only select GO terms with 4 to 200 genes from the lightnet
## mapped to them
#go_term_sizes <- apply(raw_go_matrix, 2, sum)
#go_matrix <- raw_go_matrix[, go_term_sizes >= 4 & go_term_sizes <= 200]

# Match up the gene set matrix with the prediction/query matrices:
common_genes = intersect(rownames(gene_set_matrix), colnames(best_prediction_mat))
gene_set_matrix = gene_set_matrix[common_genes, , drop = FALSE]
best_prediction_mat <- best_prediction_mat[, common_genes, drop = FALSE]
best_query_mat <- best_query_mat[, common_genes, drop = FALSE]
best_query_name_mat = best_query_name_mat[, common_genes, drop = FALSE]

# Remove GO terms from the gene set matrix without any annotations (must only be
# annotated from genes not in the set that is predicted against in the gene-level
# target prediction step).
nonzero_degree_gene_set_inds = colSums(gene_set_matrix) > 0
gene_set_matrix = gene_set_matrix[, nonzero_degree_gene_set_inds, drop = FALSE]
print(sprintf('Removed %s gene sets with no annotations', sum(!nonzero_degree_gene_set_inds)))


# Get p values and zscores for each drug --> go combination
# This will be done on a per-GO term and a per-drug basis, and
# the per-GO computations can be performed using as many
# different control types of samples as possible.  From this
# function call, one gets back all p values and z scores for
# ALL different computation methods.

# First, create a vector that splits the matrix up into
# treatments and controls
# treat_control_vec <- vector('character', dim(tp_mat)[1])
# treat_control_vec[] <- 'Treatment'
# treat_control_vec[grepl('DMSO', dimnames(tp_mat)[[1]], fixed = TRUE)] <- 'DMSO'
# treat_control_vec[grepl('Rand-by-strain', dimnames(tp_mat)[[1]], fixed = TRUE)] <- 'Rand-by-strain'


# test_set <- c(which(treat_control_vec == 'Treatment')[1:20],
#               which(treat_control_vec == 'DMSO')[1:20],
#               which(treat_control_vec == 'Rand-by-strain')[1:20]
#               )

# Select a subset for debugging purposes
# Selects 100 of each sample type (or fewer, if fewer than
# 100 of a sample type
test_inds <- do.call(c, lapply(unique(all_control_vec), function(x) {
                               inds <- which(all_control_vec == x)
                               inds[1:(min(length(inds), 100))]
}))


# Is this a test to make sure everything works, or is this the real deal?
# Look for the command line option "--test_run"
if (opt$test_run) {
    condition_inds <- test_inds
} else {
    condition_inds <- TRUE
}

# Settle the per-condition randomization question
num_rand = config_params$Required_arguments$`num_per-condition_randomizations`

# Set the seed (required input from the user)
seed = as.numeric(config_params$Required_arguments$`per-condition_randomization_seed`)
print(seed)
if (is.numeric(seed)){
    set.seed(seed)
} else if (seed == 'rand') {
    seed
} else {
    stop('specified per-condition_randomization_seed is neither numeric nor "rand"')
}

# stop('time to figure out what is going terribly wrong')

# Now, get all versions of p values and zscores!
system.time(all_pvals_zscores <- compute_zscores_pvals_by_go_and_drug(target_prediction_mat = best_prediction_mat[condition_inds,],
                                                                      go_term_mat = gene_set_matrix,
                                                                      sample_type_split_vec = all_control_vec[condition_inds],
                                                                      control_types = controls,
                                                                      types_to_predict_for = sample_types_to_predict_for,
                                                                      num_by_drug_rand = num_rand))

print(str(all_pvals_zscores))
print(controls)

print(which(is.na(all_pvals_zscores[['per_gene_set']][[controls[1]]][['pval']]), arr.ind = TRUE))
print(sum(is.na(all_pvals_zscores[['per_gene_set']][[controls[1]]][['pval']])))

print(which(is.na(all_pvals_zscores[['per_gene_set']][[controls[1]]][['zscore']]), arr.ind = TRUE))
print(sum(is.na(all_pvals_zscores[['per_gene_set']][[controls[1]]][['zscore']])))

bad_GO_term_inds = unique(which(is.na(all_pvals_zscores[['per_gene_set']][[controls[1]]][['zscore']]), arr.ind = TRUE)[, 2])
print(str(gene_set_matrix[, bad_GO_term_inds]))
print(colSums(gene_set_matrix[, bad_GO_term_inds]))

#print(which(is.na(all_pvals_zscores[['per_gene_set']][[controls[2]]][['pval']])))
#print(which(is.na(all_pvals_zscores[['per_gene_set']][[controls[2]]][['zscore']])))
# get worst per-go pvals/zscores
if (length(all_pvals_zscores[['per_gene_set']]) == 2) {
    per_gene_set_worst_case_mats <- get_worst_case_pval_zscore(all_pvals_zscores[['per_gene_set']][[controls[1]]][['pval']],
                                                               all_pvals_zscores[['per_gene_set']][[controls[1]]][['zscore']],
                                                               all_pvals_zscores[['per_gene_set']][[controls[2]]][['pval']],
                                                               all_pvals_zscores[['per_gene_set']][[controls[2]]][['zscore']],
                                                               controls[1],
                                                               controls[2]
                                                               )
} else if (length(all_pvals_zscores[['per_gene_set']]) == 1) {
    print(str(all_pvals_zscores))
    per_gene_set_worst_case_mats = list(worst_pval = all_pvals_zscores[['per_gene_set']][[controls[1]]][['pval']],
                                        worst_zscore = all_pvals_zscores[['per_gene_set']][[controls[1]]][['zscore']],
                                        control_name = matrix(controls[1], nrow = dim(all_pvals_zscores[['per_gene_set']][[controls[1]]])[1], ncol = dim(all_pvals_zscores[['per_gene_set']][[controls[1]]])[2], dimnames = dimnames(all_pvals_zscores[['per_gene_set']][[controls[1]]]))
                                        )
} else {
    stop(sprintf('Not sure how this happened, as there are %s control types when there should be either 1 or 2', length(controls)))
}

# get worst overall pvals/zscores (worst between
# per-gene-set and per-condition
overall_worst_case_mats <- get_worst_case_pval_zscore(per_gene_set_worst_case_mats$worst_pval,
                                                      per_gene_set_worst_case_mats$worst_zscore,
                                                      all_pvals_zscores$per_condition$pval,
                                                      all_pvals_zscores$per_condition$zscore,
                                                      'per_gene_set',
                                                      'per_condition'
                                                      ) 
# Specify which per-gene-set scheme was used if the worst pval/zscore
# was generated using the per-gene-set scheme
per_gene_set_indices <- overall_worst_case_mats$control_name == 'per_gene_set'
overall_worst_case_mats$control_name[per_gene_set_indices] <- per_gene_set_worst_case_mats$control_name[per_gene_set_indices]

# I need to subset the target prediction matrix so I can get the correct
# process prediction drivers back!
best_prediction_subset_mat <- best_prediction_mat[condition_inds, ][all_control_vec[condition_inds] %in% sample_types_to_predict_for, ]
best_query_subset_mat = best_query_mat[condition_inds, ][all_control_vec[condition_inds] %in% sample_types_to_predict_for, ]
best_query_name_subset_mat = best_query_name_mat[condition_inds, ][all_control_vec[condition_inds] %in% sample_types_to_predict_for, ]

# Here I will get a table of drivers of my go process predictions
gene_set_drivers_dt <- get_go_drivers(best_prediction_subset_mat, best_query_subset_mat, best_query_name_subset_mat, gene_set_matrix, cutoff = driver_cutoff)
# If the driver gene names were not specified, remove that column from the target prediction table!
if (is.null(gene_name_column)) {
    gene_set_drivers_dt[, driver_name := NULL]
}

# print(gene_set_drivers_dt[condition == 'SHANGHAI-1511_000088'])

# Set the key on the go_drivers_dt so that it can be joined with the
# process prediction data.tables
setkeyv(gene_set_drivers_dt, c('condition', 'gene_set'))

## Melt per-gene-set only predictions into a data table
#per_gene_set_worst_case_df <- lapply(per_gene_set_worst_case_mats, as.vector)
#per_gene_set_worst_case_df <- as.data.frame(per_gene_set_worst_case_df)
#
#per_gene_set_GENE_SET <- as.character(col(per_gene_set_worst_case_mats$worst_pval, as.factor = TRUE))
#per_gene_set_CONDITION <- as.character(row(per_gene_set_worst_case_mats$worst_pval, as.factor = TRUE))
#
#per_go_worst_case_dt <- data.table(condition = per_gene_set_CONDITION, gene_set = per_gene_set_GENE_SET, per_gene_set_worst_case_df)

# Melt overall only predictions into a data table
overall_worst_case_df <- lapply(overall_worst_case_mats, as.vector)
overall_worst_case_df <- as.data.frame(overall_worst_case_df)

overall_GENE_SET <- as.character(col(overall_worst_case_mats$worst_pval, as.factor = TRUE))
overall_CONDITION <- as.character(row(overall_worst_case_mats$worst_pval, as.factor = TRUE))

overall_worst_case_dt <- data.table(condition = overall_CONDITION, gene_set = overall_GENE_SET, overall_worst_case_df)


# Join the process score tables with the process driver table
#setkeyv(per_go_worst_case_dt, c('condition', 'gene_set'))
setkeyv(overall_worst_case_dt, c('condition', 'gene_set'))

#per_gene_set_worst_case_dt <- gene_set_drivers_dt[per_gene_set_worst_case_dt]
overall_worst_case_dt <- gene_set_drivers_dt[overall_worst_case_dt]

# Add a column indicating control/treatment category for each condition.
overall_worst_case_dt[, sample_type := all_control_vec[condition]]

# If the user provided a table mapping gene set IDs to interpretable gene set names,
# then add those to the final table!
if (!is.null(gene_set_name_tab)) {
#    setkey(per_gene_set_worst_case_dt, gene_set)
    setkey(overall_worst_case_dt, gene_set)
    setkey(gene_set_name_tab, gene_set)
#    per_gene_set_worst_case_dt = per_gene_set_worst_case_dt[gene_set_name_tab]
    overall_worst_case_dt = overall_worst_case_dt[gene_set_name_tab, nomatch = 0]
}

# If the user provided a table mapping condition IDs to to interpretable condition names,
# then add those to the final table!
if (!is.null(condition_name_tab)) {
#    setkey(per_gene_set_worst_case_dt, gene_set)
    setkey(overall_worst_case_dt, condition)
    setkey(condition_name_tab, condition)
#    per_gene_set_worst_case_dt = per_gene_set_worst_case_dt[gene_set_name_tab]
    overall_worst_case_dt = overall_worst_case_dt[condition_name_tab, nomatch = 0]
}


## Check size distribution of these go terms
#go_terms <- per_go_worst_case_dt[, unique(GO)]
#go_termsizes <- colSums(go_matrix[, go_terms])
#hist(go_termsizes, breaks = 200)
## Comment: I think the size distribution checks out just fine :)
## However, these terms may be slightly different from the ones
## I have predicted in the past.  I am more confident in these
## propagations, however

##### Print out some diagnostic plots describing how all of the pvals/zscores
##### relate to each other across the types of control conditions I used to
##### derive them
####
##### Need to fix speed, or DO NOT RUN UNLESS YOU HAVE EXTRA TIME!!!
##### Also, plots probably look fugly at this point
#####with(all_pvals_zscores$per_go, print_diagnostic_plots_pval_zscore_2_schemes(DMSO$pval,
#####                                                                            DMSO$zscore,
#####                                                                            `Rand-by-strain`$pval,
#####                                                                            `Rand-by-strain`$zscore,
#####                                                                            'DMSO',
#####                                                                            'Rand-by-strain',
#####                                                                            file.path(prefix, sprintf('%s_process-prediction_qc-plots', Sys.Date()))
#####                                                                            ))
#####
#####print_diagnostic_plots_pval_zscore_2_schemes(per_go_worst_case_mats$worst_pval,
#####                                             per_go_worst_case_mats$worst_zscore,
#####                                             overall_worst_case_mats$worst_pval,
#####                                             overall_worst_case_mats$worst_zscore,
#####                                             'per-GO',
#####                                             'per-drug',
#####                                             file.path(prefix, sprintf('%s_process-prediction_qc-plots', Sys.Date()))
#####                                             )
#####
#####
#
#
## Sort each prediction table by condition, pval, and zscore
#setkeyv(per_go_worst_case_dt, c('drug'))
#per_go_worst_case_dt <- per_go_worst_case_dt[, .SD[order(-worst_pval, worst_zscore, decreasing = TRUE)], by = list(drug)]
#
setkey(overall_worst_case_dt, condition)
# overall_worst_case_dt <- overall_worst_case_dt[, .SD[order(-worst_pval, worst_zscore, decreasing = TRUE)], by = list(drug)]
overall_worst_case_dt <- overall_worst_case_dt[order(worst_pval, -worst_zscore)]

# Set the order of the columns!!! (This varies depending on which other tables were provided and joined to the predictions)
print(names(overall_worst_case_dt))
cols = c('condition', 'condition_name', 'gene_set', 'gene_set_name', 'p_value', 'z_score', 'driver_query', 'driver_name', 'driver_score', 'sample_type', 'control_name')
setnames(overall_worst_case_dt, c('worst_pval', 'worst_zscore'), c('p_value', 'z_score'))
if (is.null(gene_set_name_tab)) cols = c(cols[1:(which(cols == 'gene_set_name') - 1)], cols[(which(cols == 'gene_set_name') + 1):length(cols)])
if (is.null(condition_name_tab)) cols = c(cols[1:(which(cols == 'condition_name') - 1)], cols[(which(cols == 'condition_name') + 1):length(cols)])
if (is.null(gene_name_column)) cols = c(cols[1:(which(cols == 'driver_name') - 1)], cols[(which(cols == 'driver_name') + 1):length(cols)])
setcolorder(overall_worst_case_dt, cols)

# Write out to (compressed) table-formatted file
#outfile_per_go <- gzfile(file.path(prefix, sprintf('%s_per-go_pval-zscore.txt.gz', Sys.Date())), open = 'wt')
#write.table(per_go_worst_case_dt, file = outfile_per_go, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
#close(outfile_per_go)


outfile_overall <- gzfile(output_table, open = 'wb')
write.table(overall_worst_case_dt, file = outfile_overall, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
close(outfile_overall)

# Write out a top table and a globular table for human consumption

# First, make filenames for the top table and glob top table outfiles
#out_top_per_go <- file.path(prefix, sprintf('%s_per-go_pval-zscore_top-ten.txt.gz', Sys.Date()))
#out_top_glob_per_go <- file.path(prefix, sprintf('%s_per-go_pval-zscore_top-ten_glob.txt.gz', Sys.Date()))
## Now export prediction tables for just per_go predictions
#top_ten_tab_per_go <- export_top_table(dat = per_go_worst_case_dt, outfile = out_top_per_go, select_expr = rank(worst_pval, ties.method = 'min') <= 10, split_by_expr = list(drug), order_within_by_expr = order(-worst_pval, worst_zscore, decreasing = TRUE), connection_FUN = 'gzfile')
#table_to_glob(dat = top_ten_tab_per_go, outfile = out_top_glob_per_go, by_vec = c('drug'), connection_FUN = 'gzfile')

# More filenames
#out_top_overall <- file.path(prefix, sprintf('%s_overall_pval-zscore_top-ten.txt.gz', Sys.Date()))
#out_top_glob_overall <- file.path(prefix, sprintf('%s_overall_pval-zscore_top-ten_glob.txt.gz', Sys.Date()))
## And, export the overall pval/zscore top table and glob top table files
#top_ten_tab_overall <- export_top_table(dat = overall_worst_case_dt, outfile = out_top_overall, select_expr = rank(worst_pval, ties.method = 'min') <= 10, split_by_expr = list(drug), order_within_by_expr = order(-worst_pval, worst_zscore, decreasing = TRUE), connection_FUN = 'gzfile')
#table_to_glob(dat = top_ten_tab_overall, outfile = out_top_glob_overall, by_vec = c('drug'), connection_FUN = 'gzfile')




