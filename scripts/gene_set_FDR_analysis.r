#!/usr/bin/env Rscript

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

library(optparse)
library(reshape2)
library(data.table)

# Source in libraries for gene set target prediction
TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/FDR_analysis.r'))
source(file.path(TARGET_PATH, 'lib/top_table.r'))
source(file.path(TARGET_PATH, 'lib/table_to_globular.r'))
source(file.path(TARGET_PATH, 'lib/pos_args.r'))
source(file.path(TARGET_PATH, 'lib/scott_themes.r'))

positional_arguments_list = c(gene_set_predictions = 'Gzipped file containing the gene-set-level target predictions for each condition. Must have one column each for p-values, z-scores, and the type of sample (treatment, expt_control, rand-by-strain).',
                              output_folder = 'Folder into which the results are written. Results include the gene_set_predictions file with a false discovery rate column added, as well as graphs illustrating the FDR estimation process.'
                              )

option_list = list(
    make_option('--sample_table', help = 'File with a table that maps conditions to condition-specific information'),
    make_option('--include_column', help = 'Name of True/False column in the sample_table that indicates which samples should actually be used to estimate the false discovery rate. For how we do our screens, this should be a column that is FALSE for the positive controls (we run a TON of them) and TRUE for everything else. For example, I typically make a column called "not_pos_control".'),
    make_option(c('-v', '--verbosity'), type = 'integer', default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1, is default.')
    )

final_usage_positional = get_usage_and_positional(positional_arguments_list)

parser = OptionParser(usage = final_usage_positional, option_list = option_list)
arguments = parse_args(parser, positional_arguments = length(positional_arguments_list))
opt = arguments$options
arg = arguments$args

print(opt)
print(arg)

################# Load in the required files, deal with arguments
# Read in gene_set_predictions file (can let fread determine column classes, as the
# condition ID is a concatenation of screen_name and expt_id.
gene_set_prediction_file = normalizePath(arg[1])
gene_set_prediction_tab = fread(sprintf('gzip -dc %s', gene_set_prediction_file))

# Take care of sample table argument
true_false_vec = c('True' = TRUE, 'False' = FALSE, 'TRUE' = TRUE, 'FALSE' = FALSE)
if (!is.null(opt[['sample_table']])) {
    sample_table = fread(opt[['sample_table']], colClasses = 'character')

    if (!is.null(opt[['include_column']])) {
        include_tab = sample_table[, list(condition = sprintf('%s_%s', screen_name, expt_id), include = true_false_vec[sample_table[[opt[['include_column']]]]])]
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
        message('--sample_table was specified, but --include_columns was not. Ignoring --sample_table option')
    }
} else {
    gene_set_prediction_tab_split_out = NULL
}



# Create output directory
out_dir = normalizePath(arg[2])
dir.create(out_dir, recursive = TRUE)

# Get top prediction for each condition
top_pval_zscore_dt <- get_top_prediction_per_condition(gene_set_prediction_tab)
if (!is.null(gene_set_prediction_tab_split_out)) {
    top_pval_zscore_dt_split_out <- get_top_prediction_per_condition(gene_set_prediction_tab_split_out)
}

# Add scaled counts per type of sample (Treatment, DMSO, Rand-by-strain)
top_pval_zscore_dt <- add_top_discovery_counts(top_pval_zscore_dt)

print(top_pval_zscore_dt)

# Plot discoveries vs p value
plot_discoveries_vs_pval(top_pval_zscore_dt, 'discoveries-vs-pval', out_dir)

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
plot_fdr_vs_discoveries(top_pval_zscore_dt, controls, 'fdr_vs_discoveries', out_dir)

# Split the gene set target predictions into treatment, expt_control, and random tables
gene_set_prediction_tab_treatment = gene_set_prediction_tab[sample_type == 'treatment']
gene_set_prediction_tab_expt_control = gene_set_prediction_tab[sample_type == 'expt_control']
gene_set_prediction_tab_rand_by_strain = gene_set_prediction_tab[sample_type == 'rand-by-strain']

# Write out process prediction and top process prediction files
# with FDRs attached
# And a simpler version for viewing
write_pval_zscore_dt(gene_set_prediction_tab_treatment, 'gene_set_target_prediction_treatment_with_fdr', out_dir)
write_pval_zscore_dt(gene_set_prediction_tab_expt_control, 'gene_set_target_prediction_expt_control_with_fdr', out_dir)
write_pval_zscore_dt(gene_set_prediction_tab_rand_by_strain, 'gene_set_target_prediction_rand-by-strain_with_fdr', out_dir)
if (!is.null(gene_set_prediction_tab_split_out)) {
    write_pval_zscore_dt(gene_set_prediction_tab_split_out, 'gene_set_target_prediction_split_out_with_fdr', out_dir)
}

#write_top_pval_zscore_dt(top_pval_zscore_dt, 'top_gene_set_prediction_per_condition_with_fdr', out_dir)
#if (!is.null(gene_set_prediction_tab_split_out)) {
#    write_top_pval_zscore_dt(top_pval_zscore_dt_split_out, 'top_gene_set_prediction_per_condition_split_out_with_fdr', out_dir)
#}
#write_top_pval_zscore_dt_glob_by_sample_type(top_pval_zscore_dt, 'top_gene_set_prediction_per_condition_with_fdr_glob', out_dir)
#if (!is.null(gene_set_prediction_tab_split_out)) {
#    write_top_pval_zscore_dt_glob_by_sample_type(top_pval_zscore_dt_split_out, 'top_gene_set_prediction_per_condition_split_out_with_fdr_glob', out_dir)
#}
