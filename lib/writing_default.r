#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

ext_switch <- function(con_name) {

    switch(con_name,
           file = 'txt',
           gzfile = 'txt.gz'
           )
}

get_filename_def <- function(con_name, filename_meat, prefix) {
    file.path(prefix, sprintf('%s.%s',filename_meat, ext_switch(con_name)))
}

write_tab_def <- function(tab, con_name, filename_meat, prefix, ...) {
    
    con <- match.fun(con_name)

    filename <- get_filename_def(con_name, filename_meat, prefix)
    out_con <- con(filename, open = 'wt', ...)
    write.table(tab, out_con, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
    close(out_con)

}
