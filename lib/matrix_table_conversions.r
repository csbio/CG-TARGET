#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

library(reshape2)

matrix_to_table = function(mat, row_id_tab, col_id_tab) {

    row_id_tab_with_key = copy(row_id_tab)
    row_id_tab_with_key[, row_key := sprintf('%s_%s', Strain_ID, Barcode)]
    
    col_id_tab_with_key = copy(col_id_tab)
    col_id_tab_with_key[, col_key := sprintf('%s_%s', screen_name, expt_id)]

    dimnames(mat)[1] = row_id_tab_with_key[['key']]
    dimnames(mat)[2] = col_id_tab_with_key[['key']]

    tab = melt(mat, varnames = c('row_key', 'col_key'), value.name = 'score')

    # Join the row names to the table
    setkey(tab, row_key)
    setkey(row_id_tab_with_key, key)
    tab = tab[row_id_tab_with_key]
    
    # Join the column names to the table
    setkey(tab, col_key)
    setkey(col_id_tab_with_key, key)
    tab = tab[col_id_tab_with_key]

    return(tab[, list(Strain_ID, Barcode, screen_name, expt_id, value)])
}

table_to_matrix = function(tab) {
    
    tab[, row_key := sprintf('%s_%s', Strain_ID, Barcode)]
    tab[, col_key := sprintf('%s_%s', screen_name, Strain_ID)]
    
    row_id_tab = unique(tab[, list(row_key, Strain_ID, Barcode)])
    setkey(row_id_tab, row_key)
    col_id_tab = unique(tab[, list(col_key, screen_name, expt_id)])
    setkey(col_id_tab, col_key)

    mat = acast(tab, row_key ~ col_key, value.var = 'score')
    
    tab[, row_key := NULL]
    tab[, col_key := NULL]


    row_id_tab_final = row_id_tab[dimnames(mat)[1]]
    col_id_tab_final = col_id_tab[dimnames(mat)[2]]
    dimnames(mat) = NULL

    return(list(row_id_tab = row_id_tab_final,
                col_id_tab = col_id_tab_final,
                matrix = mat
                ))
}

