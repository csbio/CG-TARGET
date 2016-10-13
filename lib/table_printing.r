#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

library(data.table)

prettify_table = function(x) {

    # Set some parameters
    min_description_width = 40
    max_total_width = 120

    x_mat = as.matrix(x)

    #print(x_mat)

    x_titles = names(x)
    
    #print(x_titles)

    # Get the column widths. This takes into account one border
    # for each cell, as well as the space padding on both sides.
    # It also takes the title of each column into account.
    col_widths_row_one = nchar(x_titles)
    col_widths_init = apply(x, 2, function(y) max(nchar(y)))
    col_widths = pmax(col_widths_row_one, col_widths_init) + 3

    #print(col_widths)

    # Determine the amount of padding (3 spaces per cell, plus
    # one extra
    #total_padding = 3 * length(col_widths) + 1

    # Get the width of the description column, which depends on
    # the total allowed width and the width of the other
    # columns.
    description_width = max(min_description_width + 3,
                            ((max_total_width - 1) - sum(col_widths[x_titles != 'description'])))

    #print(col_widths[x_titles != 'description'])
    #print(sum(col_widths[x_titles != 'description']))

    # Integrate description width back into the column widths!
    col_widths[x_titles == 'description'] = description_width

    #print(description_width)

    wrapped_description = strwrap(x_mat[, x_titles == 'description'], width = description_width - 3, simplify = FALSE)
    description_heights = vapply(wrapped_description, length, numeric(1))

    #print(wrapped_description)

    num_rows = sum(description_heights) + length(description_heights) + 3
    num_cols = length(col_widths)

    final_mat = matrix(character(), nrow = num_rows, ncol = num_cols)

    # Iterate over the original row structure, going nested into the
    # heights of each row. I therefore need a row counter for the
    # final matrix
    for (col_ind in 1:num_cols) {

        #message(1)

        # Turn the column into a list so I can easily index into it
        # and get back NA values when I need an extra row (due to
        # the wrapped description column that takes up multiple
        # rows per cell.
        #print(x_mat[, col_ind])
        wrapped_cells = strwrap(x_mat[, col_ind], width = col_widths[col_ind] - 3, simplify = FALSE)

        #message(2)

        #if (x_titles[col_ind] == 'description') {
        #    col_list = wrapped_description
        #} else {
        #    col_list = as.list(x_mat[, col_ind])
        #}
       
        # Set up the format string for all entries in the column!
        if (col_ind == num_cols) {
            format_string = sprintf('|%% %ds |', (col_widths[col_ind] - 2))
            title_string = sprintf('|%s|', paste(rep('=', col_widths[col_ind] - 1), collapse = ''))
            bottom_border_string = sprintf('|%s|', paste(rep('-', col_widths[col_ind] - 1), collapse = ''))
        } else {
            format_string = sprintf('|%% %ds ', (col_widths[col_ind] - 2))
            title_string = sprintf('|%s', paste(rep('=', col_widths[col_ind] - 1), collapse = ''))
            bottom_border_string = sprintf('|%s', paste(rep('-', col_widths[col_ind] - 1), collapse = ''))
        }

        #message(3)
        
        # First, get the title set up
        
        final_mat[1, col_ind] = title_string
        final_mat[2, col_ind] = sprintf(format_string, x_titles[col_ind])
        final_mat[3, col_ind] = title_string
    
        row_ind = 4

        # Now, loop over the rows
        for (i in seq_along(description_heights)) {

            for (j in 1:description_heights[i]) {
                
                raw_string = wrapped_cells[[i]][j]
                if (is.na(raw_string)) {
                    raw_string = ''
                }

                final_string = sprintf(format_string, raw_string)
                final_mat[row_ind, col_ind] = final_string
            
                row_ind = row_ind + 1
            }

            # At the end of each original row, add a bottom border of hyphens
            final_mat[row_ind, col_ind] = bottom_border_string
            row_ind = row_ind + 1

        }

    }
   
    #return(final_mat)

    # Now that I have the final matrix, smash it together and return the text string!
    final_vec = apply(final_mat, 1, paste, collapse = '')
    final_string = paste(final_vec, collapse = '\n')

    return(final_string)

}
