#!/usr/bin/env Rscript

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# In this script I will predict targets for the LDA-normalized dataset
# I will use Raamesh's 'dotcosine' method that normalizes the SGA matrix but not the chemgen data
# Also, I will only use the 'LightNet' genes from SGA, which have higher quality profiles (~1500)
# These profiles have degree greater than a certain threshold

library(optparse)
library(reshape2)
library(data.table)
library(yaml)

# Source in the similarity metric(s) I will use 
# NOTE: I am not computing p values
TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/cosine.r'))
source(file.path(TARGET_PATH, 'lib/dot_cosine.r'))
source(file.path(TARGET_PATH, 'lib/pos_args.r'))
source(file.path(TARGET_PATH, 'lib/filenames.r'))
source(file.path(TARGET_PATH, 'lib/datasets.r'))
source(file.path(TARGET_PATH, 'lib/dummy_dataset.r'))

positional_arguments_list = c(CONFIG_FILE = 'yaml-formatted gene-set prediction configuration file.')

option_list = list(make_option('--rand', action = 'store_true', default = FALSE, help = 'Use this flag if you want to perform the predictions for the resampled profiles')
                   )
final_usage_positional = get_usage_and_positional(positional_arguments_list)
parser = OptionParser(usage = final_usage_positional, option_list = option_list)
arguments = parse_args(parser, positional_arguments = length(positional_arguments_list))
print(arguments)
arg = arguments$args
opt = arguments$options

print(arg)

config_f = file(arg[1], 'rt')
config_params = yaml.load_file(config_f)
close(config_f)

# Read in the CG data
output_folder = config_params$Required_arguments$output_folder
if (opt$rand) {

    # Again, if the user has specified that no resampled profiles are to be generated...
    # Stop everything!!!
    zero_resampled_profiles = (config_params$Required_arguments$`per-array_resampling_scheme` == 0) | (config_params$Required_arguments$`num_per-array_resampled_profiles`== 0)
    if (zero_resampled_profiles) {
        warning('Config file specifies that no resampled profiles should be generated.\nSince no resampled profiles have been or will be generated, this\nscript to predict their gene targets using a genetic interaction\nnetwork is therefore unnecessary. Exiting now.')
        quit('no')
    }

    
    folder = get_resampled_profile_folder(output_folder)
    cg_filename = get_resampled_profile_filename(folder)
} else {
    cg_filename = config_params$Required_arguments$cg_data_table
}

cg_tab = fread(sprintf('gzip -dc %s', cg_filename), colClasses = c('character', 'character', 'character', 'character', 'numeric'))
cg_row_tab = fread(config_params$Required_arguments$cg_row_info_table, colClasses = 'character')
cg_col_tab = fread(config_params$Required_arguments$cg_col_info_table, colClasses = 'character')

# If a dummy dataset was specified, and it has experimental controls,
# then add these into the dataset!

bool_vec = c(`True` = TRUE, `False` = FALSE, `TRUE` = TRUE, `FALSE` = FALSE)
dummy_name = config_params$Options$dummy_dataset$name
if (!is.null(dummy_name)) {
    dummy_config_f = file(get_dummy_config_filename(dummy_name), 'rt')
    dummy_config_params = yaml.load_file(dummy_config_f)
    close(dummy_config_f)
    dummy_dt = fread(sprintf('gzip -dc %s', file.path(get_dummy_folder(dummy_name), dummy_config_params$cg_data_table)), colClasses = c('character', 'character','character','character','numeric'))
    dummy_col_tab = fread(file.path(get_dummy_folder(dummy_name), dummy_config_params$cg_col_info_tab), header = TRUE, colClasses = 'character')

    # If a column specifying negative experimental controls was specified, then use it
	# to filter the dummy dataset. otherwise, do not use the dummy dataset here!
    if (!is.null(config_params$Options$dummy_dataset$negative_control_column)) {
        print(unique(dummy_dt[, list(Strain_ID, Barcode)], by = NULL))
        print(unique(dummy_dt[, list(screen_name, expt_id)], by = NULL))
        message(sprintf('dummy matrix dimensions before filtering for "cols_to_include": (%s, %s)',
                      dim(unique(dummy_dt[, list(Strain_ID, Barcode)], by = NULL))[1],
                      dim(unique(dummy_dt[, list(screen_name, expt_id)], by = NULL))[1]))
        select_rows_dummy = bool_vec[dummy_col_tab[[config_params$Options$dummy_dataset$negative_control_column]]]
        dummy_col_tab = dummy_col_tab[select_rows_dummy]
    }

    # Filter the dummy interactions to just the controls
    setkeyv(dummy_col_tab, c('screen_name', 'expt_id'))
    setkeyv(dummy_dt, c('screen_name', 'expt_id'))
    dummy_col_key = dummy_col_tab[, list(screen_name, expt_id)]
    dummy_dt = dummy_dt[dummy_col_key, nomatch = 0]
    message(sprintf('Filtered dummy matrix dimensions: (%s, %s)',
                  dim(unique(dummy_dt[, list(Strain_ID, Barcode)], by = NULL))[1],
                  dim(unique(dummy_dt[, list(screen_name, expt_id)], by = NULL))[1]))

    # Create a final sample table with the dummy controls added in.
	cg_col_tab = add_dummy_controls_to_cg_sample_tab(cg_col_tab,
                                                     config_params$Options$gene_set_target_prediction$negative_control_column,
                                                     config_params$Options$gene_set_target_prediction$condition_name_column,
                                                     dummy_col_tab,
                                                     config_params$Options$dummy_dataset$negative_control_column,
                                                     config_params$Options$dummy_dataset$condition_name_column)

    # Now, combine the cg dataset with the dummy dataset controls, limiting to strains
    # found in both datasets.
    cg_strains = unique(cg_tab[, list(Barcode, Strain_ID)])
    dummy_strains = unique(dummy_dt[, list(Barcode, Strain_ID)])
    setkeyv(cg_strains, c('Barcode', 'Strain_ID'))
    setkeyv(dummy_strains, c('Barcode', 'Strain_ID'))
    strain_Barcode_intersect = cg_strains[dummy_strains, nomatch = 0]

    setkeyv(cg_tab, c('Barcode', 'Strain_ID'))
    setkeyv(dummy_dt, c('Barcode', 'Strain_ID'))
    cg_tab = rbind(cg_tab[strain_Barcode_intersect], dummy_dt[strain_Barcode_intersect])
}

message(sprintf('Final CG matrix dimensions: (%s, %s)',
                dim(unique(cg_tab[, list(Strain_ID, Barcode)], by = NULL))[1],
                dim(unique(cg_tab[, list(screen_name, expt_id)], by = NULL))[1]))

# read in the gi data
gi_info = get_gi_info(config_params$Required_arguments$gi_dataset_name, TARGET_PATH)
gi_tab = fread(sprintf('gzip -dc %s', gi_info$gi_tab))
gi_array_tab = fread(gi_info$gi_array_tab)
gi_to_cg_match_col = gi_info$array_sys_name_col

# shape the cg data into a matrix with the rownames coming from
# the column that matches the gi data match column
setkeyv(cg_tab, c('Strain_ID', 'Barcode'))
setkeyv(cg_row_tab, c('Strain_ID', 'Barcode'))
cg_tab[, condition_key := sprintf('%s_%s', screen_name, expt_id)]
cg_tab = cg_tab[cg_row_tab[, c('Strain_ID', 'Barcode', config_params$Required_arguments$cg_to_gi_match_col), with = FALSE], nomatch = 0, allow.cartesian = TRUE]
print(cg_tab)
cg_form = as.formula(sprintf('%s ~ %s', config_params$Required_arguments$cg_to_gi_match_col, 'condition_key'))
cg_mat = acast(data = cg_tab, formula = cg_form, fill = 0, fun.aggregate = mean, na.rm = TRUE, value.var = 'score')
print(str(cg_mat))

# shape the gi data into a matrix with the rownames coming from
# the column that matches the cg data match column
setkey(gi_tab, array_key)
setkey(gi_array_tab, array_key)
gi_tab = gi_tab[gi_array_tab[, c('array_key', gi_to_cg_match_col), with = FALSE], nomatch = 0, allow.cartesian = TRUE]
print(gi_tab)
gi_form = as.formula(sprintf('%s ~ %s', gi_to_cg_match_col, 'query_key'))
# while the queries (columns) in the matrix should be unique, the arrays,
# since they are defined only by a systematic gene name, could be duplicated
# and will therefore be averaged if that is the case.
gi_mat = acast(data = gi_tab, formula = gi_form, fill = 0, fun.aggregate = mean, na.rm = TRUE, value.var = 'score')
print(str(gi_mat))

# filter the matrices for genes (rows/arrays) that do not exist
# in both, and unify the row order
intersect_genes = intersect(rownames(cg_mat), rownames(gi_mat))
cg_mat <- cg_mat[intersect_genes, ]
gi_mat <- gi_mat[intersect_genes, ]

# filter out any columns in either matrix that do not have any interactions.
# this would occur after intersecting the genes common to the cg and gi
# profiles, if for a particular profile, no interactions were measured with
# the set of common genes.
# this causes problems with the similarity calculations, as each column is
# divided by its norm (divide by zero problem).
cg_mat = cg_mat[, colSums(abs(cg_mat), na.rm = TRUE) != 0]
gi_mat = gi_mat[, colSums(abs(gi_mat), na.rm = TRUE) != 0]

# predict targets (calculate similarities between columns of gi_mat and cg_mat)
# this uses the either 'cosine' similarity or 'dotcosine' similarity, the latter
# of which normalizes the columns of the gi_matrix but not the cg_matrix
predictions_mat = switch(config_params$Required_arguments$cg_gi_similarity_measure,
                         `cg-norm-dot` = dot_cos_sim(cg_mat, gi_mat, na.rm = TRUE),
                         cosine = cos_sim(cg_mat, gi_mat, na.rm = TRUE)
                         )

# melt the predictions back into a table and add the cg screen_name and expt_id
# columns back to the table!
predictions_tab = data.table(melt(predictions_mat, varnames = c('condition_key', 'query_key'), value.name = 'score'))
print(predictions_tab)
setkey(predictions_tab, condition_key)
condition_tab = unique(cg_tab[, list(condition_key, screen_name, expt_id)], by = NULL)
print(condition_tab)
setkey(condition_tab, condition_key)
predictions_tab = predictions_tab[condition_tab, nomatch = 0][, list(screen_name, expt_id, query_key, score)]
print(predictions_tab)

# open a gzipped outfile connection
out_dir = get_gene_target_folder(output_folder)
dir.create(out_dir, recursive = TRUE)
if (opt$rand) {
    out_filename = get_gene_target_prediction_resampled_filename(
                                out_dir,
                                config_params$Required_arguments$`per-array_resampling_scheme`,
                                config_params$Required_arguments$`num_per-array_resampled_profiles`,
                                config_params$Required_arguments$`per-array_resampling_seed`,
                                config_params$Required_arguments$cg_gi_similarity_measure
                                )
} else {
    out_filename = get_gene_target_prediction_filename(out_dir,
                                                       config_params$Required_arguments$cg_gi_similarity_measure
                                                       )
}

pred_out <- gzfile(out_filename, open = 'wb')

# Write out to file
write.table(predictions_tab, pred_out, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
close(pred_out)

