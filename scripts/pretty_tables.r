#!/usr/bin/env Rscript

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

library(optparse)
library(reshape2)
library(data.table)
library(tools)

# Source in libraries for gene set target prediction
TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/pos_args.r'))

positional_arguments_list = c(table_file = 'Gzipped file containing the gene-set-level target predictions for each condition. Must be sorted by ascending p-value, followed by descending z-score.'
                              )

option_list = list(
    make_option(c('-v', '--verbosity'), type = 'integer', default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1, is default.')
    )

final_usage_positional = get_usage_and_positional(positional_arguments_list)

parser = OptionParser(usage = final_usage_positional, option_list = option_list)
arguments = parse_args(parser, positional_arguments = length(positional_arguments_list))
opt = arguments$options
arg = arguments$args

print(opt)
print(arg)

# Load in the file to be *prettied up*!!!
infile = arg[1]

# Get filenames together!
folder = file_path_sans_ext(file_path_sans_ext(infile))
dir.create(folder)
f_all_pred_grouped_by_cond = file.path(folder, 'condition-grouped_top-pred-sorted.txt.gz')
f_top_ten_pred_grouped_by_cond = file.path(folder, 'condition-grouped_top-ten-pred-sorted.txt.gz')
f_top_pred_per_cond = file.path(folder, 'top-pred-per-condition.txt.gz')
f_all_pred_grouped_by_cond_glob = file.path(folder, 'condition-grouped_top-pred-sorted_glob.txt.gz')
f_top_ten_pred_grouped_by_cond_glob = file.path(folder, 'condition-grouped_top-ten-pred-sorted_glob.txt.gz')

tab = fread(sprintf('zcat %s', infile))
colorder = names(tab)

# Get the top prediction for each condition
# Since the predictions are already ordered by p-value and z-score, I don't
# have to reorder them!
top_pred_tab_sorted = tab[, .SD[1], by = list(condition)][order(p_value, -z_score)]
top_cond_sorted = top_pred_tab_sorted[['condition']]

# Prepare the massive table to be reordered
# (Remember, this will not change the relative order of rows
# with the same condition - very useful!)
tab = setkey(tab, condition)

# Here, ordering the table so that all predictions for the condition with the best
# prediction are first, followed by all predictions for the condition with the
# second best prediction, etc.
tab_ordered_by_top_pred = tab[top_cond_sorted]

# Also for prettiness, grab the top ten predictions for each condition
tab_ordered_by_top_pred_top_ten = tab_ordered_by_top_pred[, .SD[1:min(c(10, nrow(.SD)))], by = condition]

# Reset the column order of these new tables to match that of the original table
setcolorder(top_pred_tab_sorted, colorder)
setcolorder(tab_ordered_by_top_pred, colorder)
setcolorder(tab_ordered_by_top_pred_top_ten, colorder)

# Write the "long form" tables to files
all_pred_grouped_by_cond_con = gzfile(f_all_pred_grouped_by_cond, 'wt')
top_ten_pred_grouped_by_cond_con = gzfile(f_top_ten_pred_grouped_by_cond, 'wt')
top_pred_per_cond_con = gzfile(f_top_pred_per_cond, 'wt')

write.table(tab_ordered_by_top_pred, all_pred_grouped_by_cond_con, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(tab_ordered_by_top_pred_top_ten, top_ten_pred_grouped_by_cond_con, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(top_pred_tab_sorted, top_pred_per_cond_con, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

close(all_pred_grouped_by_cond_con)
close(top_ten_pred_grouped_by_cond_con )
close(top_pred_per_cond_con)

# Get the columns for the globbed output tables
glob_cols = colorder[!colorder %in% c('condition', 'condition_name')]

# Write out the globbed files!
# For all predictions
all_pred_glob_con = gzfile(f_all_pred_grouped_by_cond_glob, 'wt')
setkey(tab_ordered_by_top_pred, condition)
for (cond in top_cond_sorted) {
    write(sprintf('condition: %s', cond), file = all_pred_glob_con)
    if ('condition_name' %in% names(tab_ordered_by_top_pred)) {
        write(sprintf('condition_name: %s', tab_ordered_by_top_pred[cond][['condition_name']][1]), file = all_pred_glob_con)
    }
    write('---', file = all_pred_glob_con)
    write.table(tab_ordered_by_top_pred[cond][, glob_cols, with = FALSE], file = all_pred_glob_con, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    write('\n\n-------------------------\n\n', file = all_pred_glob_con)
}
close(all_pred_glob_con)

# For top ten predictions
top_ten_glob_con = gzfile(f_top_ten_pred_grouped_by_cond_glob, 'wt')
setkey(tab_ordered_by_top_pred_top_ten, condition)
for (cond in top_cond_sorted) {
    write(sprintf('condition: %s', cond), file = top_ten_glob_con)
    if ('condition_name' %in% names(tab_ordered_by_top_pred_top_ten)){
        write(sprintf('condition_name: %s', tab_ordered_by_top_pred_top_ten[cond][['condition_name']][1]), file = top_ten_glob_con)
    }
    write('---', file = top_ten_glob_con)
    write.table(tab_ordered_by_top_pred_top_ten[cond][, glob_cols, with = FALSE], file = top_ten_glob_con, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    write('\n\n-------------------------\n\n', file = top_ten_glob_con)
}
close(top_ten_glob_con)

