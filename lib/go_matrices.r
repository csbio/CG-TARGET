#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

library(GO.db)
library(graph)
library(RBGL)
library(org.Sc.sgd.db)
library(foreach)
library(reshape2)
library(data.table)


get_gene_go_bp_is_a_propagated_matrix <- function(orf_universe) {
    
    get_go_2_go_map <- function() {
        
        # Open a connection to the GO database
        go <- GO_dbconn()

        # Create a temp table of the go bp parent relationships,
        # but ONLY the is_a relationships!
        create_is_a_bp_parent_table <-
            'CREATE TEMP TABLE go_bp_parents_is_a AS SELECT _id, _parent_id FROM go_bp_parents WHERE relationship_type = "is_a";'
        dbSendQuery(go, create_is_a_bp_parent_table)

        #     dbGetQuery(go, 'SELECT _id from go_bp_parents_is_a LIMIT 10;')

        # Get the go term names from the 'children' or 'from' column
        go_children <- dbGetQuery(go, 'SELECT gt.go_id
                        FROM go_term AS gt,
                             go_bp_parents_is_a AS gbp
                             WHERE gt._id = gbp._id
                             ;')$go_id
        
                             #     dbGetQuery(go, 'SELECT _parent_id from go_bp_parents_is_a LIMIT 10;')

        # Get the go term names from the 'parent' or 'to' column
        go_parents <- dbGetQuery(go, 'SELECT gt.go_id
                        FROM go_term AS gt,
                             go_bp_parents_is_a AS gbp
                             WHERE gt._id = gbp._parent_id
                             ;')$go_id

        # Drop the temporary is_a table
        dbSendQuery(go, 'DROP TABLE go_bp_parents_is_a;')

        # Get list-style mapping and return
        tapply(go_parents, INDEX = go_children, FUN = function(x) x)
        

    }

    get_orf_2_go_map <- function() {
        # Get a table of gene to GO annotations for S. cerevisiae,
        # for the orf_universe!
        orf_universe_go_df <- select(org.Sc.sgd.db, keys = orf_universe, columns = 'GO', keytype = 'ORF')
        orf_universe_go_bp_df <- orf_universe_go_df[orf_universe_go_df$ONTOLOGY == 'BP' & !is.na(orf_universe_go_df$GO), c('ORF', 'GO')]

        # Get list-style mapping and return
        tapply(orf_universe_go_bp_df$GO, INDEX = orf_universe_go_bp_df$ORF, FUN = function(x) x)
        
    }

    prop_par <- function(x, go_2_go) {

        if (any(x %in% names(go_2_go))) {
#            message(paste('passed:', paste(x, collapse = ', ')))
#            message('\n')

            parents <- unlist(go_2_go[x], use.names = FALSE)
            parents_filtered <- parents[parents != 'all']
            return(unique(c(x,
                            prop_par(parents_filtered, go_2_go)
                            )))
        } else {
#            message(paste('failed:', paste(x, collapse = ', ')))
#            message('\n')
            return(character(0))

        }
    }
    
    
    propagate_orf_2_go_map <- function(orf_2_go, go_2_go) {

        orf_2_go_prop <- lapply(orf_2_go, prop_par, go_2_go)
        names(orf_2_go_prop) <- names(orf_2_go)

        orf_2_go_prop

    }

    map_list_2_df <- function(map_list) {
        
        map_dfs <- lapply(names(map_list), function(x) {
                          data.frame(ORF = x, GO = map_list[[x]])}
        )

        as.data.frame(rbindlist(map_dfs))

    }


    orf_2_go_map_2_matrix <- function(orf_2_go) {

        orf_2_go_df <- map_list_2_df(orf_2_go)
#        print(head(orf_2_go_df))
        names(orf_2_go_df) <- c('ORF', 'GO')
        orf_2_go_df$edge_exists <- 1
        acast(orf_2_go_df, formula = ORF ~ GO, value.var = 'edge_exists', fill = 0)

    }


    # Here we actually execute all of the functions!
    orf_2_go_map <- get_orf_2_go_map()
    go_2_go_map <- get_go_2_go_map()
    
    orf_2_go_prop_map <- propagate_orf_2_go_map(orf_2_go_map, go_2_go_map)
    
    orf_2_go_map_2_matrix(orf_2_go_prop_map)

}

get_only_go_bp_is_a_propagated_matrix <- function(orf_universe) {

    get_go_2_go_map <- function() {
        
        # Open a connection to the GO database
        go <- GO_dbconn()

        # Create a temp table of the go bp parent relationships,
        # but ONLY the is_a relationships!
        create_is_a_bp_parent_table <-
            'CREATE TEMP TABLE go_bp_parents_is_a AS SELECT _id, _parent_id FROM go_bp_parents WHERE relationship_type = "is_a";'
        dbSendQuery(go, create_is_a_bp_parent_table)

        #     dbGetQuery(go, 'SELECT _id from go_bp_parents_is_a LIMIT 10;')

        # Get the go term names from the 'children' or 'from' column
        go_children <- dbGetQuery(go, 'SELECT gt.go_id
                        FROM go_term AS gt,
                             go_bp_parents_is_a AS gbp
                             WHERE gt._id = gbp._id
                             ;')$go_id
        
                             #     dbGetQuery(go, 'SELECT _parent_id from go_bp_parents_is_a LIMIT 10;')

        # Get the go term names from the 'parent' or 'to' column
        go_parents <- dbGetQuery(go, 'SELECT gt.go_id
                        FROM go_term AS gt,
                             go_bp_parents_is_a AS gbp
                             WHERE gt._id = gbp._parent_id
                             ;')$go_id

        # Drop the temporary is_a table when the function exits
        on.exit(dbSendQuery(go, 'DROP TABLE go_bp_parents_is_a;'))

        # Get list-style mapping and return
        tapply(go_parents, INDEX = go_children, FUN = function(x) x)
        

    }

    get_orf_2_go_map <- function() {
        # Get a table of gene to GO annotations for S. cerevisiae,
        # for the orf_universe!
        orf_universe_go_df <- select(org.Sc.sgd.db, keys = orf_universe, columns = 'GO', keytype = 'ORF')
        orf_universe_go_bp_df <- orf_universe_go_df[orf_universe_go_df$ONTOLOGY == 'BP' & !is.na(orf_universe_go_df$GO), c('ORF', 'GO')]

        # Get list-style mapping and return
        tapply(orf_universe_go_bp_df$GO, INDEX = orf_universe_go_bp_df$ORF, FUN = function(x) x)
        
    }

    prop_par <- function(x, go_2_go) {

        if (any(x %in% names(go_2_go))) {
#            message(paste('passed:', paste(x, collapse = ', ')))
#            message('\n')

            parents <- unlist(go_2_go[x], use.names = FALSE)
            parents_filtered <- parents[parents != 'all']
            return(unique(c(x,
                            prop_par(parents_filtered, go_2_go)
                            )))
        } else {
#            message(paste('failed:', paste(x, collapse = ', ')))
#            message('\n')
            return(character(0))

        }
    }
    
    get_go_init <- function(orf_2_go) {

        unique(unlist(orf_2_go, use.names = FALSE))

    }

    propagate_go_2_go_map <- function(go_init, go_2_go) {

        go_2_go_prop <- lapply(go_init, prop_par, go_2_go)
        names(go_2_go_prop) <- go_init 

        go_2_go_prop

    }

    map_list_2_df <- function(map_list) {
        
        map_dfs <- lapply(names(map_list), function(x) {
                          data.frame(ORF = x, GO = map_list[[x]])}
        )

        as.data.frame(rbindlist(map_dfs))

    }


    orf_2_go_map_2_matrix <- function(orf_2_go) {

        orf_2_go_df <- map_list_2_df(orf_2_go)
#        print(head(orf_2_go_df))
        names(orf_2_go_df) <- c('ORF', 'GO')
        orf_2_go_df$edge_exists <- 1
        acast(orf_2_go_df, formula = ORF ~ GO, value.var = 'edge_exists', fill = 0)

    }


    # Here we actually execute all of the functions!
    orf_2_go_map <- get_orf_2_go_map()
    go_2_go_map <- get_go_2_go_map()
   
    go_children <- get_go_init(orf_2_go_map)
    go_2_go_prop_map <- propagate_go_2_go_map(go_children, go_2_go_map)
    
    orf_2_go_map_2_matrix(go_2_go_prop_map)

}


go_orf_overlap_hyper_p <- function(go1, go2, orf_2_go_mat) {

    x <- orf_2_go_mat

    #     message(sprintf('go1: %s', go1))
    #     message(sprintf('go2: %s', go2))

    num_orfs_drawn_in_go1 <- sum(intersect(x[, go1], x[, go2]))
    num_orfs_in_go1 <- sum(x[, go1])
    num_orfs_not_in_go1 <- sum(abs(1 - x[, go1]))
    num_orfs_drawn <- sum(x[, go2])

    # I have to do 1 - p for consistency with other
    # overlap metrics
    1 - phyper(q = num_orfs_drawn_in_go1 - 1,
           m = num_orfs_in_go1,
           n = num_orfs_not_in_go1,
           k = num_orfs_drawn,
           lower.tail = FALSE
           )

}


go_orf_overlap_jaccard <- function(go1, go2, orf_2_go_mat) {

    x <- orf_2_go_mat
    #     message(sprintf('go1: %s', go1))
    #     message(sprintf('go2: %s', go2))
    
    go1_bool_map <- as.logical(x[, go1])
    go2_bool_map <- as.logical(x[, go2])

    # Return Jaccard similarity index
    sum(go1_bool_map & go2_bool_map) / sum(go1_bool_map | go2_bool_map)

}


go_orf_max_overlap <- function(go1, go2, orf_2_go_mat) {

    x <- orf_2_go_mat
    #     message(sprintf('go1: %s', go1))
    #     message(sprintf('go2: %s', go2))
    
    go1_bool_map <- as.logical(x[, go1])
    go2_bool_map <- as.logical(x[, go2])

    # Return maximum fraction of overlap between the two go terms
    sum(go1_bool_map & go2_bool_map) / min(sum(go1_bool_map), sum(go2_bool_map))

}
