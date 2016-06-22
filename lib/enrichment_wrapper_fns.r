# I will place at least one wrapper function into this script
# It will wrap the functions I have written to perform
# superfast enrichment tests

get_enrichments <- function(results_table, condition, variables, values, test_interesting, mapping_table, map_from, map_to) {

    # results_table: the table that contains condition, variables, and values columns
    # condition: name of the column that groups together results on which enrichment is tested (example, each different sample in a microarray experiment could be a condition)
    # variables: name of the column that specifies which variables were measured in the experiment (example, the gene names in a microarray experiment).  These variables will be mapped to GO terms or other lists for enrichment testing
    # values: name of the column that possesses the measured values for each of the variables under each condition
    # test_interesting: an expression that tests every element in 'value' to see if it is in the 'interesting' set
    # mapping_table: the table that maps the variables to meaningful lists, such as GO terms
    # map_from: the name of the column in the mapping table whose values match those in the 'variables' column
    # map_to: the name of the column in the mapping table that contains informative list names
    
    # Enrichment is tested between each pair of the values found in 'condition' and 'map_to'
    
    # Load in external libraries
    require(data.table)
    source('../lib/obo_tools.r')

    # Convenience function to pass character strings into data.table
    z <- function(x) {
        quote(parse(text = x))
    }

    # Substitute all column name and expression variables!
    cond <- substitute(condition)
    vars <- substitute(variables)
    vals <- substitute(values)
    test <- substitute(test_interesting)
    from <- substitute(map_from)
    to <- substitute(map_to)

    # Ensure that the tables are data.tables and that they do not modify
    # existing data.tables in the user's environment
    results_table <- copy(as.data.table(results_table))
    mapping_table <- copy(as.data.table(mapping_table))

    # Set keys on the results and mapping tables
    setkeyv(results_table, deparse(cond))
    setkeyv(mapping_table, deparse(from))

    #     print(mapping_table)

    # Get the 'variable universe.'  Basically, just the names of all of the variables
    # measured in this experiment.
    var_universe <- unique(results_table[[deparse(vars)]])
    var_universe_size <- length(var_universe)
    #     print(mapping_table[, eval(from)])

    # Ensure the mapping is unique (otherwise false enrichments can occur if duplicate
    # genes map to the same go term)
    # Also, ensure that only the variables reported on in the study
    # are used for mapping to the informative lists.
    # If there are extra variables that map to the informative lists,
    # then the statistics get messed up
    # (example: if there are 5 genes in an informative list but only three were
    # investigated in a study, the informative list should only include those three genes.
    mapping_table <- unique(mapping_table, by = c(deparse(from), deparse(to)))[var_universe, , nomatch = 0]

    # Get the number of items that map to each meaningful list
    # Only use the variables that were reported on to start the mapping
    # Otherwise the terms are too large and the math is not correct
    #     mapping_table[, num_vars_in_term := length(eval(from))), by = list(eval(to))]
    eval(bquote(mapping_table[, num_vars_in_term := .N, by = .(to)]))

    #     print(mapping_table)
    
         #     print(mapping_table[, list(num_vars_in_term = .N), by = list(eval(to))])

    # Debugging
    #     print(results_table[eval(test)])

    #     print(results_table[eval(test)][J(unique(eval(cond))), .SD[, c(.SD, list(num_int_vars = length(eval(vars))))]])

    # Perform the enrichment test!  (at compiled C speed!!!!!)
    ### VERY IMPORTANT NOTE!!!!!!!!
    ### When performing the p value calculation, it is performed using the overlap - 1
    # and any other parameters that depend on the overlap have 1 added to them accordingly.
    # This has to do with how the p value is calculated and prevents overestimation of significance.
    # For example, if the overlap is 6, the p value is calculated as the probability that
    # and overlap GREATER THAN 6 would occur.  However, we want to take the probability of
    # our observation into account as well, so we would set the overlap to 5 instead of 6.
    # Also, no enrichments will be reported if the overlap is only one variable - this will never
    # be a relevant overlap!
    enrich_expr <- bquote(results_table[J(unique(.(cond))), .SD[.(test), c(.SD, list(num_int_vars = .N))][, .SD, keyby = .(vars)][mapping_table, , nomatch = 0][, c(.SD, list(num_int_vars_in_term = .N)), keyby = .(to)][num_int_vars_in_term > 1, c(.SD, list(overlapping = paste(.(vars), collapse = '|'))), by = .(to)][, unique(.SD)], by = .EACHI][, .SD, keyby = list(.(cond), .(to))][, list(num_int_vars_in_term, num_vars_in_term, num_vars_not_in_term = (var_universe_size - num_vars_in_term), num_int_vars, p_enrich = phyper(num_int_vars_in_term - 1, num_vars_in_term, (var_universe_size - num_vars_in_term), num_int_vars, lower.tail = FALSE), overlapping = overlapping), by = list(.(cond),.(to))])
    #     print(enrich_expr)

    eval(enrich_expr)

}

