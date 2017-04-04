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
cg_dt <- fread(sprintf('gzip -dc %s', config_params$Required_arguments$cg_data_table), header = TRUE, colClasses = c('character', 'character','character','character','numeric'))
cg_row_tab = fread(config_params$Required_arguments$cg_row_info_table, header = TRUE, colClasses = 'character')
cg_col_tab = fread(config_params$Required_arguments$cg_col_info_table, header = TRUE, colClasses = 'character')
print(cg_dt)

# Get the info tabs ready to roll
setkeyv(cg_row_tab, c('Strain_ID', 'Barcode'))
setkeyv(cg_col_tab, c('screen_name', 'expt_id'))

# If a cols_to_include column was specified, filter the row info table, which
# is then used to filter the data themselves
if (!is.null(config_params$Options$`per-array_randomization`$`per-array_randomization_include_column`)) {
    #I don't think I need this flag anymore...
    #if (arg$`per-array_randomization_include_column` == 'TRUE') {
    #    stop('The cols_to_include flag is present, but no argument is specified!')
    #}
    print(unique(cg_dt[, list(Strain_ID, Barcode)], by = NULL))
    print(unique(cg_dt[, list(screen_name, expt_id)], by = NULL))
    message(sprintf('Matrix dimensions before filtering for "cols_to_include": (%s, %s)',
                  dim(unique(cg_dt[, list(Strain_ID, Barcode)], by = NULL))[1],
                  dim(unique(cg_dt[, list(screen_name, expt_id)], by = NULL))[1]))
    bool_vec = c(`True` = TRUE, `False` = FALSE, `TRUE` = TRUE, `FALSE` = FALSE)
    select_rows = bool_vec[cg_col_tab[[config_params$Options$`per-array_randomization`$`per-array_randomization_include_column`]]]
    cg_col_tab = cg_col_tab[select_rows]
}
# Here is where the data are filtered. The filtering happens regardless of the
# presence of the cols_to_include flag.
setkeyv(cg_col_tab, c('screen_name', 'expt_id'))
setkeyv(cg_dt, c('screen_name', 'expt_id'))
cg_col_key = cg_col_tab[, list(screen_name, expt_id)]
cg_dt = cg_dt[cg_col_key, nomatch = 0]
message(sprintf('Final matrix dimensions: (%s, %s)',
              dim(unique(cg_dt[, list(Strain_ID, Barcode)], by = NULL))[1],
              dim(unique(cg_dt[, list(screen_name, expt_id)], by = NULL))[1]))

setkeyv(cg_dt, c('Strain_ID', 'Barcode'))

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
    rand_dt <- cg_dt[, list(screen_name = combined_screen_name, expt_id = sprintf('%06d', 1:n), score = sample(score, n, replace = TRUE)), by = list(Strain_ID, Barcode)]
} else if (config_params$Required_arguments$`per-array_resampling_scheme` == 2) {
    rand_dt <- cg_dt[, list(screen_name = combined_screen_name, expt_id = sprintf('%06d', 1:n), score = rnorm(n, mean(score, na.rm = TRUE), sd = sd(score, na.rm = TRUE))), by = list(Strain_ID, Barcode)]
} else {
    stop('"per-array_resampling_scheme" was given, but value was not "0", "1", or "2"!')
}
setcolorder(rand_dt, c('Strain_ID', 'Barcode', 'screen_name', 'expt_id', 'score'))

out_dir = get_resampled_profile_folder(config_params$Required_arguments$output_folder)
dir.create(out_dir, recursive = TRUE)
out_filename = get_resampled_profile_filename(out_dir)
of = gzfile(out_filename, 'wb')
write.table(rand_dt, of, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
close(of)


