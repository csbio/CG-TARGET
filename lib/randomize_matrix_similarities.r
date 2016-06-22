# This script will utilize the similarity functions in this folder (and perhaps other folders)
# and compute empirical p values for the results they generate

# The functions build up incrementally and are based on vectors.  This is still efficient when
# working with matrics because you can perform your fast matrix-based similarity calculations
# on matrices where each column is a randomization of one vector.

# The implementation concerned with calculating empirical p values for similarities between all
# columns to two matrices is parallelized using the 'parallel package'

# Move this comment someplace else...
# Obtain an empirical p value by shuffling the row labels on matrix 1,
# calculating similarity, and comparing to the distribution of randomized p values.

# First, load in some required libraries
library(foreach)
library(iterators)
library(parallel)
library(doParallel)

# Define a function to combine the results from each eigendrug into two pval matrices
# (pos and neg).  Make sure it takes in more than 2 inputs!
combine_pvals <- function(...) {
    
    # Takes in as many inputs as necessary
    # Each input value should be a list with 'pos' and 'neg' vectors
    # I slice the pos and neg out of each one and put them in their own lists
    pvals <- list(...)
    list(pos = do.call(rbind, lapply(pvals, `[[`, 'pos')), neg = do.call(rbind, lapply(pvals, `[[`, 'neg')))
}

# Start a progress bar (mostly for debugging purposes...)
# pb <- txtProgressBar(title = "Predicting eigendrug targets", style = 3)
# n <- length(comps_to_analyze)


# Define a function that returns p values given the following:
# 1) A vector to be randomized; 2) a function to compute results of the randomizations;
# and 3) a vector of observed results to compare to the randomized results

get_empirical_pvals <- function(obs_data_vec, obs_res_vec, n, FUN, ...) {
    
    compute_randomized_results <- match.fun(FUN)

    get_randomizations <- function() {
        # Returns a matrix in which each column is a randomization of obs_data_vec
        replicate(n, sample(obs_data_vec, replace = FALSE), simplify = TRUE)

    }

    get_randomized_results <- function(randomizations) {
        # FUN should be able to take in a matrix of randomizations (each column
        # is a randomization) and compute the relevant result for each of those
        # randomized columns/vectors.  

        compute_randomized_results(randomizations, ...)
        
    }

    get_randomized_pvals <- function(randomized_results) {
        
        # The following assumes that the result of calling 'get_randomized_results' is organized
        # such that each randomization is now a row.  If this is not the default of your
        # function, FUN will need to wrap that function and generate the correct data structure.

        # Sweep out a logical matrix indicating which randomized results were more extreme,
        # either positively or negatively, than the observed results (obs_rec_vec)
        #    message('getting p values part 1')
        bool_p_pos <- sweep(randomized_results, 2, as.vector(obs_res_vec), '>=')
        bool_p_neg <- sweep(randomized_results, 2, as.vector(obs_res_vec), '<=')

        # Get the sum of each column of the pos and neg boolean cutoff matrices
        #    message('getting p values part 2')
        sum_p_pos <- apply(bool_p_pos, 2, sum)
        sum_p_neg <- apply(bool_p_neg, 2, sum)

        # Compute p values by dividing by the total number of randomizations (N)
        #         message('getting p values part 3')
        p_pos <- sum_p_pos / n
        p_neg <- sum_p_neg / n

        return(list(pos = p_pos, neg = p_neg))
    }

    randomizations <- get_randomizations()
    randomized_results <- get_randomized_results(randomizations)
    randomized_pvals <- get_randomized_pvals(randomized_results)

    #     message('randomizations')
    #     message(str(randomizations))
    #     message('randomized_results')
    #     message(str(randomized_results))
    #     message('randomized_pvals')
    #     message(str(randomized_pvals))

    return(randomized_pvals)
}


# Define a function that generates empirical p values in the case of calculating similarities
# between all columns of two matrices.  It is parallelized
get_colwise_emp_pvals <- function(data_mat, res_mat, FUN, n, ...) {
    
    # Register parallel backend
    # Just a performance note: I achieve ~3x faster performance when using 6 cores
    # Obviously these cores are all on the same machine
    registerDoParallel(cores = 6)

    # For the moment, '...' contains any remaining arguments to the inner functions that I can't remember
    # I can't remember if 'n' and 'FUN' will be passed in without error...

    # data_mat is the original data that will have something done to it (FUN) in order to generate
    # the corresponding ROW in res_mat.

    # Set up iterators over columns of data and results matrices
    data_cols <- iter(data_mat, by = 'column')
    res_cols <- iter(res_mat, by = 'row')

    emprical_pvals <- foreach(data_col = data_cols, res_col = res_cols, i = icount(), .combine = combine_pvals, .multicombine = TRUE) %dopar% {

        data_col <- as.vector(data_col)
        res_col <- as.vector(res_col)

        message(paste('Computing p values for column', i, 'predictions', sep = ' '))

        get_empirical_pvals(data_col, res_col, n, FUN, ...)

    }

    # And give them the correct names again
    dimnames(emprical_pvals$pos)[[1]] <- dimnames(res_mat)[[1]]
    dimnames(emprical_pvals$neg)[[1]] <- dimnames(res_mat)[[1]]

    return(emprical_pvals)
}



