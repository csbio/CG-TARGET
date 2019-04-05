#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################


# In this script, I write two functions that take in
# a vector of test statistics or other values and a
# vector containing a distribution against which to
# test the values in the first vector. The first
# function resembles what is already in my 
# chemical genomics process prediction pipeline,
# and the second function is the "new" way of
# doing it that should be faster and use less
# memory.

# This is basically the old function, although it was not
# nicely wrapped like this :)
#two_vector_sig_old = function(test_vec, distrib_vec) {
#    num_controls = length(distrib_vec)
#    vapply(test_vec, function(x) {
#           sum(x <= distrib_vec) / num_controls
#           }, numeric(1))
#}

empirical_pval_factory = function(alternative = c('greater', 'less')) {

    if (alternative == 'greater') {
        compare_fn = `>=`
    } else if (alternative == 'less') {
        compare_fn = `<=`
    } else {
        stop('The argument for "alternative" must be in c("greater", "less")')
    }
        

    f =  function(test_vec, distrib_vec) {

        i_test_vec = order(test_vec, decreasing = TRUE)
        i_distrib_vec = order(distrib_vec, decreasing = TRUE)

        num_distrib = length(distrib_vec)
        num_test = length(test_vec)

        test_sig = numeric(num_test)

        i = 1
        j = 1
        while (i <= num_distrib & j <= num_test) {
            # Weird way to do the comparison, but it works and makes
            # it so the "greater" and "less" functions are always the
            # same. Functional programming!!!
            if (compare_fn(distrib_vec[i_distrib_vec[i]], test_vec[i_test_vec[j]])) {
                i = i + 1        
            } else {
                test_sig[i_test_vec[j]] = i - 1
                j = j + 1
            }

        }
        
        if (j <= num_test) {
            test_sig[i_test_vec[j:num_test]] = i - 1
        }

        test_sig / num_distrib
    }

    return(f)
}

empirical_pval_greater = empirical_pval_factory('greater')
empirical_pval_less = empirical_pval_factory('less')

#empirical_pval_greater = function(test_vec, distrib_vec) {
#    i_test_vec = order(test_vec, decreasing = TRUE)
#    i_distrib_vec = order(distrib_vec, decreasing = TRUE)
#
#    num_distrib = length(distrib_vec)
#    num_test = length(test_vec)
#
#    test_sig = numeric(num_test)
#
#    i = 1
#    j = 1
#    while (i <= num_distrib & j <= num_test) {
#        if (distrib_vec[i_distrib_vec[i]] >= test_vec[i_test_vec[j]]) {
#            i = i + 1        
#        } else if (distrib_vec[i_distrib_vec[i]] < test_vec[i_test_vec[j]]) {
#            test_sig[i_test_vec[j]] = i - 1
#            j = j + 1
#        }
#
#    }
#    
#    if (j =< num_test) {
#        test_sig[i_test_vec[j:num_test]] = i - 1
#    }
#
#    test_sig / num_distrib
#
#}