#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This file contains two functions to take a table of data
# and write it back out to a new file in a globular, more
# human-readable form.

# The first function is main function that "globs" the table
# and writes it to a file, the second is a wrapper for the
# first that reads in a table and then export the globbed
# format

# After using this with great success for some time, I have
# decided to remove the first `setkeyv` call.  This allows
# me to construct the globs purely on the row order of the
# table that is passed in.

library(data.table)

table_to_glob <- function(dat, outfile, by_vec, connection_FUN = 'file', ...) {
    
    # Make sure we have a data.table!
    dat_dt <- data.table(dat)

    # Get a table of all unique keys
    entries <- unique(dat_dt[, by_vec, with = FALSE], by = by_vec) 
    
    # Open a file connection
    con_FUN <- match.fun(connection_FUN)
    out_con <- con_FUN(outfile, open = 'wt', ...)

    # Now write out to the file, split up by the desired columns
    dat_dt[, {
        for (i in seq_along(.BY)) {
            write(paste(by_vec[i], ': ', .BY[[i]], sep = ''), file = out_con)
        }
        write('---', file = out_con)
        write.table(.SD, file = out_con, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE, ...)
        write('\n\n-------------------------\n\n', file = out_con)
    }, 
    by = by_vec]

    close(out_con)
}

read_table_to_glob <- function(infile, outfile, by_vec, connection_FUN = 'file', sep = '\t', header = TRUE, quote = '', comment.char = '', as.is = TRUE, ...) {

    dat <- read.table(infile, sep = sep, header = header, quote = quote, comment.char = comment.char, as.is = as.is, ...)

    table_to_glob(dat = dat, outfile = outfile, by_vec = by_vec, connection_FUN = connection_FUN, ...)

}
