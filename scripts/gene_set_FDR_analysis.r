#!/usr/bin/env Rscript

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

library(optparse)
library(reshape2)
library(data.table)
library(yaml)

# Source in libraries for gene set target prediction
TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/FDR_analysis.r'))
source(file.path(TARGET_PATH, 'lib/top_table.r'))
source(file.path(TARGET_PATH, 'lib/table_to_globular.r'))
source(file.path(TARGET_PATH, 'lib/pos_args.r'))
source(file.path(TARGET_PATH, 'lib/scott_themes.r'))
source(file.path(TARGET_PATH, 'lib/filenames.r'))
source(file.path(TARGET_PATH, 'lib/dummy_dataset.r'))

positional_arguments_list = c(CONFIG_FILE = 'yaml-formatted gene-set prediction configuration file.')
option_list = list(make_option('--test_run', action = 'store_true', default = FALSE, help = 'Set this flag if you want to run the FDR estimation script on gene_set target predictions produced using the "--test_run" flag.'))

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

# Create output directories, which are based on the name of the
# gene sets used in the analysis, to allow different versions
# of the results to be generated against different gene sets.
output_folder = config_params$Required_arguments$output_folder
gene_set_name = config_params$Required_arguments$gene_set_name
final_tables_dir = get_final_results_tabs_folder(output_folder, gene_set_name)
final_plots_dir = get_final_results_plots_folder(output_folder, gene_set_name)
dir.create(final_tables_dir, recursive = TRUE)
dir.create(final_plots_dir, recursive = TRUE)

################# Load in the required files
# Read in gene_set_predictions file (can let fread determine column classes, as the
# condition ID is a concatenation of screen_name and expt_id.)
gene_set_prediction_folder = get_gene_set_target_prediction_folder(output_folder)
# Since I now have two versions of the gene_set_prediction_filename naming scheme,
# I look for the most recent version first, then the older version, while issuing
# a message stating that reading the old scheme does not take termsize into
# account.
gene_set_prediction_filename_v1 = get_gene_set_target_prediction_filename_v1(
                   gene_set_prediction_folder,
                   config_params$Required_arguments$`per-array_resampling_scheme`,
                   config_params$Required_arguments$`num_per-array_resampled_profiles`,
                   config_params$Required_arguments$`per-array_resampling_seed`,
                   config_params$Required_arguments$cg_gi_similarity_measure,
                   config_params$Required_arguments$`per-condition_randomization_seed`,
                   config_params$Required_arguments$`num_per-condition_randomizations`,
                   config_params$Required_arguments$gene_set_name,
                   opt$test_run
                   )
gene_set_prediction_filename_v2 = get_gene_set_target_prediction_filename_v2(
                   gene_set_prediction_folder,
                   config_params$Required_arguments$`per-array_resampling_scheme`,
                   config_params$Required_arguments$`num_per-array_resampled_profiles`,
                   config_params$Required_arguments$`per-array_resampling_seed`,
                   config_params$Required_arguments$cg_gi_similarity_measure,
                   config_params$Required_arguments$`per-condition_randomization_seed`,
                   config_params$Required_arguments$`num_per-condition_randomizations`,
                   config_params$Required_arguments$gene_set_name,
                   config_params$Required_arguments$min_term_size,
                   config_params$Required_arguments$max_term_size,
		   config_params$Required_arguments$alternative,
                   opt$test_run
                   )
if (file.exists(gene_set_prediction_filename_v2)) {
    gene_set_prediction_tab = fread(sprintf('gzip -dc %s', gene_set_prediction_filename_v2))
} else if (file.exists(gene_set_prediction_filename_v1)) {
    gene_set_prediction_tab = fread(sprintf('gzip -dc %s', gene_set_prediction_filename_v1))
    message('PLEASE NOTE: the process-level prediction file you are about to read in\nwas generated using an older version of CG-TARGET that did not specify\nthe gene set term sizes in the filename. If you generated the process-\nlevel prediction file using this version of CG-TARGET, please ensure that\n"predict_gene-set_target.r" ran successfuly.')
} else {
    stop(sprintf('The process-level prediction file at the location below does not exist.\nPlease ensure that the previous step in the pipeline, "predict_gene-set_targets.r",\nran successfully and/or the config parameters did not change between steps.\nMissing file:\n%s', gene_set_prediction_filename_v2))
}

# Take care of sample table argument
true_false_vec = c('True' = TRUE, 'False' = FALSE, 'TRUE' = TRUE, 'FALSE' = FALSE)
sample_table_filename = config_params$Required_arguments$cg_col_info_table
sample_table = fread(sample_table_filename, header = TRUE, colClasses = 'character')

# Deal with the dummy dataset sample table, as in previous scripts
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
		select_rows_dummy = true_false_vec[dummy_col_tab[[config_params$Options$dummy_dataset$negative_control_column]]]
		dummy_col_tab = dummy_col_tab[select_rows_dummy]
	}

	# Create a final sample table with the dummy controls added in.
	sample_table = add_dummy_controls_to_cg_sample_tab(sample_table,
													   config_params$Options$gene_set_target_prediction$negative_control_column,
													   config_params$Options$gene_set_target_prediction$condition_name_column,
													   dummy_col_tab,
													   config_params$Options$dummy_dataset$negative_control_column,
													   config_params$Options$dummy_dataset$condition_name_column)
}

######## Here, which samples are included in the analysis are determined
######## Also, it is determined here if there are conditions to "split out,"
######## which would be conditions that are neither negative controls or
######## deemed as "treatments" in the config file.

treatment_column = config_params$Options$results_export$treatment_column

if (!is.null(treatment_column)) {
  
    # First, check if treatment column is in the cg column table (sample table)
    if (! treatment_column %in% names(sample_table)) {
        stop(treatment_column %in% names(sample_table),
             sprintf('Options : results_export : treatment_column : "%s"\nis not in the cg_col_info_table found at:\n%s',
                     treatment_column, sample_table_filename)
             )
    }

    # First, convert the treatment_column to TRUE/FALSE if it's not already
    treatment_vec = true_false_vec[sample_table[[treatment_column]]]

    # Then grab the negative control column and also convert it to TRUE/FALSE
    # If the control column exists, then the negative controls were used in
    # gene_set target prediction and should be included **in the FDR
    # calculations**, so we do a union of treatment and negative control
    # conditions. If it does not exist, then the only conditions to include
    # **in the FDR calculations** are the treatment conditions. Any conditions
    # that are not specified as either "treatment" (using the treatment
    # column) or "negative control" (from the control column) are then
    # excluded from FDR calculations but still exported with an FDR in the
    # "split-out" output table.
    control_column = config_params$Options$gene_set_target_prediction$negative_control_column
    if (!is.null(control_column)) {
        control_vec = true_false_vec[sample_table[[control_column]]]
        include_vec = treatment_vec | control_vec
    } else {
        include_vec = treatment_vec
    }

    include_tab = sample_table[, list(condition = sprintf('%s_%s', screen_name, expt_id), include = include_vec)]
    gene_set_prediction_conditions = unique(gene_set_prediction_tab[['condition']])
    conditions_to_completely_exclude = include_tab[['condition']][!(include_tab[['condition']] %in% gene_set_prediction_conditions)]
    conditions_to_split_out = include_tab[['condition']][include_tab[['condition']] %in% gene_set_prediction_conditions & !(include_tab[['include']])]
    conditions_for_FDR = setdiff(gene_set_prediction_conditions, union(conditions_to_completely_exclude, conditions_to_split_out))
    print(length(gene_set_prediction_conditions))
    print(length(conditions_to_completely_exclude))
    print(length(conditions_to_split_out))
    print(length(conditions_for_FDR))

    gene_set_prediction_tab_split_out = gene_set_prediction_tab[condition %in% conditions_to_split_out]
    gene_set_prediction_tab = gene_set_prediction_tab[condition %in% conditions_for_FDR]

} else {
    gene_set_prediction_tab_split_out = NULL
    message('"treatment_column" was not specified. Everything that is not a negative control will be considered a "treatment."')
}


# This is specific to FDR scheme 1
fdr_scheme = config_params$Required_arguments$FDR_estimation_scheme

if (fdr_scheme == 1) {
    # Get top prediction for each condition
    top_pval_zscore_dt <- get_top_prediction_per_condition(gene_set_prediction_tab)
    if (!is.null(gene_set_prediction_tab_split_out)) {
        top_pval_zscore_dt_split_out <- get_top_prediction_per_condition(gene_set_prediction_tab_split_out)
    }

    # Add scaled counts per type of sample (Treatment, DMSO, Rand-by-strain)
    top_pval_zscore_dt <- add_top_discovery_counts(top_pval_zscore_dt)

    print(top_pval_zscore_dt)

    # Plot discoveries vs p value
    plot_discoveries_vs_pval(top_pval_zscore_dt, 'discoveries_vs_pval', final_plots_dir)

    # Calculate FDRs and get a p value to FDR map
    pval_fdr_dt <- get_pval_fdr_mapping(top_pval_zscore_dt)

    ## Join FDR to top and all process prediction tables
    sample_types <- unique(top_pval_zscore_dt[['sample_type']])
    treatment <- 'treatment'
    controls <- sample_types[sample_types != treatment]
    top_pval_zscore_dt <- add_fdr(top_pval_zscore_dt, pval_fdr_dt, controls)
    gene_set_prediction_tab <- add_fdr(gene_set_prediction_tab, pval_fdr_dt, controls)
    if (!is.null(gene_set_prediction_tab_split_out)) {
        gene_set_prediction_tab_split_out <- add_fdr(gene_set_prediction_tab_split_out, pval_fdr_dt, controls)
    }

    # Plot FDR vs discoveries per FDR type
    plot_fdr_vs_discoveries(top_pval_zscore_dt, controls, 'fdr_vs_discoveries', final_plots_dir)

} else if (fdr_scheme == 2) {

    col_names_raw = names(gene_set_prediction_tab)
    col_order = c('p_value', 'BH_FDR', col_names_raw[col_names_raw != 'p_value'])
    gene_set_prediction_tab[, BH_FDR := p.adjust(p_value, 'BH'), by = condition]
    setcolorder(gene_set_prediction_tab, col_order)
    if (!is.null(gene_set_prediction_tab_split_out)) {
        gene_set_prediction_tab_split_out[, BH_FDR := p.adjust(p_value, 'BH'), by = condition]
        setcolorder(gene_set_prediction_tab_split_out, col_order)
    }

    
} else {

    stop('Valid values for the FDR_estimation_scheme parameter are 1 or 2 -- not "%s", as specified in the given config file.', fdr_scheme)

}

# Split the gene set target predictions into treatment, expt_control, and random tables
gene_set_prediction_tab_treatment = gene_set_prediction_tab[sample_type == 'treatment']
gene_set_prediction_tab_expt_control = gene_set_prediction_tab[sample_type == 'expt_control']
gene_set_prediction_tab_rand_by_strain = gene_set_prediction_tab[sample_type == 'rand-by-strain']

# Write out process prediction and top process prediction files
# with FDRs attached
# And a simpler version for viewing
write_pval_zscore_dt(gene_set_prediction_tab_treatment, 'gene_set_target_prediction_treatment_with_fdr', final_tables_dir)
write_pval_zscore_dt(gene_set_prediction_tab_expt_control, 'gene_set_target_prediction_expt_control_with_fdr', final_tables_dir)
write_pval_zscore_dt(gene_set_prediction_tab_rand_by_strain, 'gene_set_target_prediction_rand-by-strain_with_fdr', final_tables_dir)
if (!is.null(gene_set_prediction_tab_split_out)) {
    write_pval_zscore_dt(gene_set_prediction_tab_split_out, 'gene_set_target_prediction_split_out_with_fdr', final_tables_dir)
}

#write_top_pval_zscore_dt(top_pval_zscore_dt, 'top_gene_set_prediction_per_condition_with_fdr', out_dir)
#if (!is.null(gene_set_prediction_tab_split_out)) {
#    write_top_pval_zscore_dt(top_pval_zscore_dt_split_out, 'top_gene_set_prediction_per_condition_split_out_with_fdr', out_dir)
#}
#write_top_pval_zscore_dt_glob_by_sample_type(top_pval_zscore_dt, 'top_gene_set_prediction_per_condition_with_fdr_glob', out_dir)
#if (!is.null(gene_set_prediction_tab_split_out)) {
#    write_top_pval_zscore_dt_glob_by_sample_type(top_pval_zscore_dt_split_out, 'top_gene_set_prediction_per_condition_split_out_with_fdr_glob', out_dir)
#}
