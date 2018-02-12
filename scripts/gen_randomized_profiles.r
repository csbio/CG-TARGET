#!/usr/bin/env Rscript

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# In this script I will read in a table of chemical-genetic profiles
# and create a large number (50,000) of random profiles.  For each
# strain, I sample a random chemical-genetic interaction between that
# strain and a random drug condition.  This way, I destroy all
# profile-specific strain covariance, but maintain the properties
# of each individual strain across all drug conditions.  If there
# are general patterns that occur in a large amount of conditions,
# these will come through in the randomized profiles.

library(optparse)
library(data.table)
library(yaml)

# Source in libraries for gene set target prediction
TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/pos_args.r'))
source(file.path(TARGET_PATH, 'lib/filenames.r'))

positional_arguments_list = c(CONFIG_FILE = 'yaml-formatted gene-set prediction configuration file.')
final_usage_positional = get_usage_and_positional(positional_arguments_list)
parser = OptionParser(usage = final_usage_positional, option_list = list())
arguments = parse_args(parser, positional_arguments = length(positional_arguments_list))
arg = arguments$args

print(arg)

config_f = file(arg[1], 'rt')
config_params = yaml.load_file(config_f)
close(config_f)

# Stop everything if the per-array_resampling_scheme is 0, or if the number
# of resampled profiles to generate is 0!
zero_resampled_profiles = (config_params$Required_arguments$`per-array_resampling_scheme` == 0) | (config_params$Required_arguments$`num_per-array_resampled_profiles`== 0)
if (zero_resampled_profiles) {
    warning('Config file specifies that no resampled profiles should be generated.\nThis is script is therefore unnecessary.\nExiting now.')
    quit('no')
}

# Read in as a data table
cg_dt <- fread(sprintf('gzip -dc %s', config_params$Required_arguments$cg_data_table), header = TRUE,
               colClasses = 'character')
cg_dt[, score := as.numeric(score)]
cg_row_tab = fread(config_params$Required_arguments$cg_row_info_table, header = TRUE, colClasses = 'character')
cg_col_tab = fread(config_params$Required_arguments$cg_col_info_table, header = TRUE, colClasses = 'character')
print(cg_dt)

# Determine if format is BEAN-counter version 1 or 2, and adjust accordingly
all_cols_v1 <- c('Strain_ID', 'Barcode', 'screen_name', 'expt_id', 'score')
all_cols_v2 <- c('Strain_ID', 'screen_name', 'expt_id', 'score')

if (all(all_cols_v1 %in% names(cg_dt))) {
  strain_id_cols <- c('Strain_ID', 'Barcode')
} else if (all(all_cols_v2 %in% names(cg_dt))) {
  strain_id_cols <- 'Strain_ID'
} else {
  stop(sprintf('The following columns must be present in the cg_data_table:\n"%s", and optionally "Barcode"\ncg_data_table location: %s\n', paste(all_cols_v2, collapse = '", "'),
               config_params$Required_arguments$cg_data_table), call. = FALSE)
}

# Get the info tabs ready to roll
setkeyv(cg_row_tab, strain_id_cols)
setkeyv(cg_col_tab, c('screen_name', 'expt_id'))

# If a cols_to_include column was specified, filter the col info table, which
# is then used to filter the data themselves
bool_vec = c(`True` = TRUE, `False` = FALSE, `TRUE` = TRUE, `FALSE` = FALSE)
if (!is.null(config_params$Options$`per-array_randomization`$`per-array_randomization_include_column`)) {
    print(unique(cg_dt[, strain_id_cols], by = NULL))
    print(unique(cg_dt[, list(screen_name, expt_id)], by = NULL))
    message(sprintf('Matrix dimensions before filtering for "cols_to_include": (%s, %s)',
                  dim(unique(cg_dt[, strain_id_cols], by = NULL))[1],
                  dim(unique(cg_dt[, list(screen_name, expt_id)], by = NULL))[1]))
    select_rows = bool_vec[cg_col_tab[[config_params$Options$`per-array_randomization`$`per-array_randomization_include_column`]]]
    cg_col_tab = cg_col_tab[select_rows]
}

# Here is where the data are filtered. The filtering happens regardless of the
# presence of the cols_to_include flag.
setkeyv(cg_col_tab, c('screen_name', 'expt_id'))
setkeyv(cg_dt, c('screen_name', 'expt_id'))
cg_col_key = cg_col_tab[, list(screen_name, expt_id)]
cg_dt = cg_dt[cg_col_key, nomatch = 0]
message(sprintf('Filtered matrix dimensions: (%s, %s)',
              dim(unique(cg_dt[, strain_id_cols], by = NULL))[1],
              dim(unique(cg_dt[, list(screen_name, expt_id)], by = NULL))[1]))

# If a dummy dataset was specified, load that in
dummy_name = config_params$Options$dummy_dataset$name
if (!is.null(dummy_name)) {
    dummy_config_f = file(get_dummy_config_filename(dummy_name), 'rt')
    dummy_config_params = yaml.load_file(dummy_config_f)
    close(dummy_config_f)
    dummy_dt = fread(sprintf('gzip -dc %s', file.path(get_dummy_folder(dummy_name), dummy_config_params$cg_data_table)), colClasses = c('character', 'character','character','character','numeric'))
    dummy_col_tab = fread(file.path(get_dummy_folder(dummy_name), dummy_config_params$cg_col_info_tab), header = TRUE, colClasses = 'character')

    # Filter the dummy dataset sample table!
    if (!is.null(config_params$Options$dummy_dataset$`per-array_randomization_include_column`)) {
        print(unique(dummy_dt[, strain_id_cols], by = NULL))
        print(unique(dummy_dt[, list(screen_name, expt_id)], by = NULL))
        message(sprintf('Dummy matrix dimensions before filtering for "cols_to_include": (%s, %s)',
                      dim(unique(dummy_dt[, strain_id_cols], by = NULL))[1],
                      dim(unique(dummy_dt[, list(screen_name, expt_id)], by = NULL))[1]))
        select_rows_dummy = bool_vec[dummy_col_tab[[config_params$Options$dummy_dataset$`per-array_randomization_include_column`]]]
        dummy_col_tab = dummy_col_tab[select_rows_dummy]
    }


    # Do the filtering for the dummy dataset.
    setkeyv(dummy_col_tab, c('screen_name', 'expt_id'))
    setkeyv(dummy_dt, c('screen_name', 'expt_id'))
    dummy_col_key = dummy_col_tab[, list(screen_name, expt_id)]
    dummy_dt = dummy_dt[dummy_col_key, nomatch = 0]
    message(sprintf('Filtered dummy matrix dimensions: (%s, %s)',
                  dim(unique(dummy_dt[, strain_id_cols], by = NULL))[1],
                  dim(unique(dummy_dt[, list(screen_name, expt_id)], by = NULL))[1]))

    # Here is where I combine dummy and "real" data, just for generating
    # the resampled profiles. Remember to match the strains!
    cg_strains = unique(cg_dt[, rev(strain_id_cols)])
    dummy_strains = unique(dummy_dt[, rev(strain_id_cols)])
    setkeyv(cg_strains, rev(strain_id_cols))
    setkeyv(dummy_strains, rev(strain_id_cols))
    strain_barcode_intersect = cg_strains[dummy_strains, nomatch = 0]

    setkeyv(cg_dt, rev(strain_id_cols))
    setkeyv(dummy_dt, rev(strain_id_cols))
    cg_dt = rbind(cg_dt[strain_barcode_intersect], dummy_dt[strain_barcode_intersect])

}
   
message(sprintf('Final matrix dimensions: (%s, %s)',
                dim(unique(cg_dt[, strain_id_cols], by = NULL))[1],
                dim(unique(cg_dt[, list(screen_name, expt_id)], by = NULL))[1]))

setkeyv(cg_dt, strain_id_cols)

# set seed, so this is reproducible, if there is a seed specified
# If not, then it is basically random each time!
seed = config_params$Required_arguments$`per-array_resampling_seed`
if (is.numeric(seed)){
    set.seed(seed)
} else if (seed == 'rand') {
    seed
} else {
    stop('specified per-array_resampling_seed is neither numeric nor "rand"')
}

# Randomize my data!
# To clarify, I am taking each strain and randomly sampling n
# values with replacement from it.  Each sample is assigned
# successively to rand-000001, rand-000002, etc.
n <- config_params$Required_arguments$`num_per-array_resampled_profiles`
combined_screen_name = sprintf('%s_%s', paste(unique(cg_dt[['screen_name']]), collapse = ','), 'rand')
if (config_params$Required_arguments$`per-array_resampling_scheme` == 1) {
    rand_dt <- cg_dt[, list(screen_name = combined_screen_name, expt_id = sprintf('%06d', 1:n), score = sample(score, n, replace = TRUE)), by = strain_id_cols]
} else if (config_params$Required_arguments$`per-array_resampling_scheme` == 2) {
    rand_dt <- cg_dt[, list(screen_name = combined_screen_name, expt_id = sprintf('%06d', 1:n), score = rnorm(n, mean(score, na.rm = TRUE), sd = sd(score, na.rm = TRUE))), by = strain_id_cols]
} else {
    stop('"per-array_resampling_scheme" was given, but value was not "0", "1", or "2"!')
}
setcolorder(rand_dt, c(strain_id_cols, 'screen_name', 'expt_id', 'score'))

out_dir = get_resampled_profile_folder(config_params$Required_arguments$output_folder)
dir.create(out_dir, recursive = TRUE)

out_filename = get_resampled_profile_filename(out_dir)
of = gzfile(out_filename, 'wb')
write.table(rand_dt, of, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
close(of)


