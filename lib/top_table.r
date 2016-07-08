#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This file contains two functions that take a table
# from one of my analyses (target predictions, GO enrichments)
# and exports a table containing only rows that meet a certain
# condition.
# 
# It will select rows within each group defined by `split_by_expr`
# and sort the results within that group by `order_within_by_expr`,
# which must be an expression of the form `order(...)`.
#
# export_top_table also returns the top table so it can be used
# by other R code, especially the table_to_glob function
#

library(data.table)

export_top_table <- function(dat, outfile, select_expr, split_by_expr, order_within_by_expr, connection_FUN = 'gzfile', ...) {

    # Make the selecting expressions usable by data.table
    select_ex <- substitute(select_expr)
    split_ex <- substitute(split_by_expr)
    order_ex <- substitute(order_within_by_expr)

    # Load in data file
    dat_dt <- data.table(dat)
    
    # Select only the rows that meet the criteria
    # given by `select_expr`
    dat_subset <- eval(bquote(dat_dt[, .SD[.(select_ex)][.(order_ex)], by = .(split_ex)]))

    # Open a file connection
    con_FUN <- match.fun(connection_FUN)
    out_con <- con_FUN(outfile, open = 'wt', ...)
    
    # Write the table out to either a normal text or gzipped file
    write.table(dat_subset, file = out_con, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE, ...)
    close(out_con)

    return(dat_subset)

}

read_table_export_top_table <- function(infile, outfile, select_expr, split_by_expr, order_within_by_expr, connection_FUN = 'gzfile', sep = '\t', header = TRUE, quote = '', comment.char = '', as.is = TRUE, ...) {
    
    dat <- read.table(infile, sep = sep, header = header, quote = quote, comment.char = comment.char, as.is = as.is, ...)

    export_top_table(dat = dat, outfile = outfile, select_expr = select_expr, split_by_expr = split_by_expr, order_within_by_expr = order_within_by_expr, connection_FUN = connection_FUN, ...)

}
