#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script contains a function, get_go_drivers,
# that will take in a target prediction matrix and
# a go term matrix (it will make sure the dimension
# labels match up), and export a list of the
# drivers of each conditions's predictions against
# all go terms.  Driver status is determined by
# passing a user-defined threshold.
#
# Note: it wants gene dim labels as orfs, and will
# convert to common names for convenience

library(data.table)

get_gene_set_drivers <- function(tp_mat, tp_query_mat, tp_query_name_mat, gene_set_mat, cutoff) {

    # Make sure tp_mat columns and gene_set_mat rows line up
    # (this should already be the case, but it never
    # hurts)
    common_genes <- intersect(dimnames(tp_mat)[[2]], dimnames(gene_set_mat)[[1]])

    tp_mat <- tp_mat[, common_genes]
    tp_query_mat = tp_query_mat[, common_genes]
    tp_query_name_mat = tp_query_name_mat[, common_genes]

    gene_set_mat <- gene_set_mat[common_genes, ]

    print(str(tp_mat))
    print(tp_mat[1:5, 1:10])
    print(str(tp_query_mat))
    print(tp_query_mat[1:5, 1:10])
    print(tp_query_mat[(dim(tp_query_mat)[1]-4):(dim(tp_query_mat)[1]), (dim(tp_query_mat)[2]-9):(dim(tp_query_mat)[2])])
    print(str(gene_set_mat))

    print(str(tp_query_name_mat))

    gene_set_bool_mat <- gene_set_mat
    mode(gene_set_bool_mat) <- 'logical'

    ## Get iterator over gene_set_term indices
    #i_gene_set_indices <- iter(gene_set_bool_mat, by = 'col')
    gene_sets = dimnames(gene_set_bool_mat)[[2]]
    num_gene_sets = dim(gene_set_bool_mat)[2]
    
    
    ###### For each gene set, compute each condition's driver target predictions ######
    # Create final data structure first!
    driver_tab_list = vector('list', num_gene_sets)

    # Another for loop?! Whaaaaa?!
    for (i in 1:num_gene_sets) {

        # Send a message to keep track of term driver progress
        cat(sprintf('assembling the query drivers for gene set %s/%s\r', i, num_gene_sets))
        #message(sprintf('assembling the query drivers for gene set %s/%s', i, num_gene_sets))

        # Get the indices for the current gene set
        gene_set_indices = gene_set_bool_mat[, i]
        #print(gene_set_indices)
        
        #message('0')

        # Get subset of target prediction matrix corresponding only to one go term
        tp_subset <- tp_mat[, gene_set_indices, drop = FALSE]
        #message('1')
        tp_query_subset = tp_query_mat[, gene_set_indices, drop = FALSE]
        #message('2')
        tp_query_name_subset = tp_query_name_mat[, gene_set_indices, drop = FALSE]
        #message('3')

        # Get indices of subsetted tp that pass the cutoff
        tp_subset_passes_cutoff <- tp_subset >= cutoff

        # Get vectors of condition, gene target, and corresponding target prediction
        # score (gene set prediction is already defined within the loop)
        # print(str(tp_subset))
        condition = as.character(row(tp_subset, as.factor = TRUE)[tp_subset_passes_cutoff])
        driver_query = as.character(tp_query_subset[tp_subset_passes_cutoff])
        driver_name = as.character(tp_query_name_subset[tp_subset_passes_cutoff])
        driver_score = tp_subset[tp_subset_passes_cutoff]
        #print(str(driver_query))

        
        # If no drugs have drivers which pass the cutoff, then return an empty data.table
        if (length(condition) == 0) {
            gene_set = character(0)
        } else {
            gene_set = gene_sets[i]
        }
        
        # Assemble into a data.table
        drivers_dt <- data.table(condition, gene_set, driver_query, driver_name, driver_score)
        setkeyv(drivers_dt, c('condition', 'gene_set'))

        
        # Get and return nicely formatted drivers and scores
        #driver_tab_list[[i]] = drivers_dt[, list(driver_query = paste(driver_query[order(driver_score, decreasing = TRUE)], collapse = '|'),
        #                                    driver_name = ifelse(all(is.na(driver_name)), NA_character_, paste(driver_name[order(driver_score, decreasing = TRUE)], collapse = '|')),
        #                                    driver_score = paste(driver_score[order(driver_score, decreasing = TRUE)], collapse = '|')
        #                                    ), by = list(condition, gene_set)
        #                                  ]
        driver_tab_list[[i]] = drivers_dt[, { driver_score_order = order(driver_score, decreasing = TRUE);
                                              list(driver_query = paste(driver_query[driver_score_order], collapse = '|'),
                                                   driver_name = ifelse(all(is.na(driver_name)), NA_character_, paste(driver_name[driver_score_order], collapse = '|')),
                                                   driver_score = paste(driver_score[driver_score_order], collapse = '|')
                                                   ) }, by = list(condition, gene_set)
                                          ]
        #print(driver_tab_list[[i]])
    }

    rbindlist(driver_tab_list)
}



    ## For each gene set, compute each condition's driver target predictions
    #foreach(gene_set_indices = i_gene_set_indices, gene_set = dimnames(gene_set_mat)[[2]], i = icount(), .combine = 'rbind', .maxcombine = 100000) %do% {
    #    
    #    # Send a message to keep track of term driver progress
    #    message(sprintf('assembling the query drivers for gene set %s/%s', i, num_gene_sets))

    #    # Convert gene_set_term_indices to a vector!
    #    gene_set_indices <- as.vector(gene_set_indices)
    #    
    #    message('0')

    #    # Get subset of target prediction matrix corresponding only to one go term
    #    tp_subset <- tp_mat[, gene_set_indices, drop = FALSE]
    #    message('1')
    #    tp_query_subset = tp_query_mat[, gene_set_indices, drop = FALSE]
    #    message('2')
    #    tp_query_name_subset = tp_query_name_mat[, gene_set_indices, drop = FALSE]
    #    message('3')


    #    # Get indices of subsetted tp that pass the cutoff
    #    tp_subset_passes_cutoff <- tp_subset >= cutoff

    #    # Get vectors of condition, gene target, and corresponding target prediction
    #    # score (GO prediction is already defined within the loop)
    #    # print(str(tp_subset))
    #    condition <- as.character(row(tp_subset, as.factor = TRUE)[tp_subset_passes_cutoff])
    #    driver_query <- as.character(tp_query_subset[tp_subset_passes_cutoff])
    #    driver_name = as.character(tp_query_name_subset[tp_subset_passes_cutoff])
    #    driver_score <- tp_subset[tp_subset_passes_cutoff]
    #    #print(str(driver_query))

    #    
    #    # If no drugs have drivers which pass the cutoff, then return an empty data.table
    #    if (length(condition) == 0) {

    #        gene_set <- character(0)
    #    } 

    #    # Assemble into a data.table
    #    drivers_dt <- data.table(condition, gene_set, driver_query, driver_name, driver_score)
    #    setkeyv(drivers_dt, c('condition', 'gene_set'))

    #    
    #    # Get and return nicely formatted drivers and scores
    #    drivers_dt[,list(driver_query = paste(driver_query[order(driver_score, decreasing = TRUE)], collapse = '|'),
    #                     driver_name = ifelse(all(is.na(driver_name)), NA_character_, paste(driver_name[order(driver_score, decreasing = TRUE)], collapse = '|')),
    #                     driver_score = paste(driver_score[order(driver_score, decreasing = TRUE)], collapse = '|')
    #                     ), by = list(condition, gene_set)
    #               ]
    #    #drivers_dt[, { driver_score_order = order(driver_score, decreasing = TRUE);
    #    #               list(driver_query = paste(driver_query[driver_score_order], collapse = '|'),
    #    #                    driver_name = ifelse(all(is.na(driver_name)), NA_character_, paste(driver_name[driver_score_order], collapse = '|')),
    #    #                    driver_score = paste(driver_score[driver_score_order], collapse = '|')
    #    #                    ) }, by = list(condition, gene_set)
    #    #           ]
    #}

#}

get_go_drivers <- function(tp_mat, tp_query_mat, tp_query_name_mat, go_mat, cutoff) {

    # Make sure tp_mat columns and go_mat rows line up
    # (this should already be the case, but it never
    # hurts)
    common_genes <- intersect(dimnames(tp_mat)[[2]], dimnames(go_mat)[[1]])

    tp_mat <- tp_mat[, common_genes]
    tp_query_mat = tp_query_mat[, common_genes]
    tp_query_name_mat = tp_query_name_mat[, common_genes]

    go_mat <- go_mat[common_genes, ]

    print(str(tp_mat))
    print(tp_mat[1:5, 1:10])
    print(str(tp_query_mat))
    print(tp_query_mat[1:5, 1:10])
    print(tp_query_mat[(dim(tp_query_mat)[1]-4):(dim(tp_query_mat)[1]), (dim(tp_query_mat)[2]-9):(dim(tp_query_mat)[2])])
    print(str(go_mat))

    go_bool_mat <- go_mat
    mode(go_bool_mat) <- 'logical'

    # Get iterator over go_term indices
    i_go_term_indices <- iter(go_bool_mat, by = 'col')
    n = ncol(go_mat)

    # For each go term, compute each condition's driver target predictions
    foreach(go_term_indices = i_go_term_indices, gene_set = dimnames(go_mat)[[2]], i = icount(), .combine = 'rbind', .maxcombine = 100000) %dopar% {
        
        # Send a message to keep track of term driver progress
        message(sprintf('assembling the query drivers for gene set %s/%s', i, n))

        # Convert go_term_indices to a vector!
        go_term_indices <- as.vector(go_term_indices)

        # Get subset of target prediction matrix corresponding only to one go term
        tp_subset <- tp_mat[, go_term_indices, drop = FALSE]
        tp_query_subset = tp_query_mat[, go_term_indices, drop = FALSE]
        tp_query_name_mat = tp_query_name_mat[, go_term_indices, drop = FALSE]

        # Get indices of subsetted tp that pass the cutoff
        tp_subset_passes_cutoff <- tp_subset >= cutoff

        # Get vectors of condition, gene target, and corresponding target prediction
        # score (GO prediction is already defined within the loop)
        # print(str(tp_subset))
        condition <- as.character(row(tp_subset, as.factor = TRUE)[tp_subset_passes_cutoff])
        driver_query <- as.character(tp_query_subset[tp_subset_passes_cutoff])
        driver_name = as.character(tp_query_name_mat[tp_subset_passes_cutoff])
        driver_score <- tp_subset[tp_subset_passes_cutoff]
        #print(str(driver_query))

        # If no drugs have drivers which pass the cutoff, then return an empty data.table
        if (length(condition) == 0) {

            gene_set <- character(0)
        } 

        # Assemble into a data.table
        drivers_dt <- data.table(condition, gene_set, driver_query, driver_name, driver_score)
        setkeyv(drivers_dt, c('condition', 'gene_set'))

        # Get and return nicely formatted drivers and scores
        drivers_dt[,list(driver_query = paste(driver_query[order(driver_score, decreasing = TRUE)], collapse = '|'),
                         driver_name = ifelse(all(is.na(driver_name)), NA_character_, paste(driver_name[order(driver_score, decreasing = TRUE)], collapse = '|')),
                         driver_score = paste(driver_score[order(driver_score, decreasing = TRUE)], collapse = '|')
                         ), by = list(condition, gene_set)
                   ]
        #drivers_dt[, { driver_score_order = order(driver_score, decreasing = TRUE);
        #               list(driver_query = paste(driver_query[driver_score_order], collapse = '|'),
        #                    driver_name = ifelse(all(is.na(driver_name)), NA_character_, paste(driver_name[driver_score_order], collapse = '|')),
        #                    driver_score = paste(driver_score[driver_score_order], collapse = '|')
        #                    ) }, by = list(condition, gene_set)
        #           ]
    }

}
