#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

library(ggplot2)
library(data.table)

# Source in libraries specific to this part of the script!
TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/empirical_pval.R'))

# target_prediction_mat should have drugs/conditions as rows and target prediction scores as columns
# go_term_mat should have the GO terms as columns and the predicted target genes as rows (and 1's in
#   the matrix if the target gene is in a go term.

# One note: I have recently changed all of the .maxcombine arguments to 100000,
# in an attempt to speed up combining all of my results once the calculations
# have completed.
best_score_per_col_group = function(mat, group_vec) {

    groups = unique(group_vec)
    num_conds = nrow(mat)
    final_mat = matrix(numeric(), nrow = num_conds, ncol = length(groups), dimnames = list(rownames(mat), groups))
    which_mat = matrix(character(), nrow = num_conds, ncol = length(groups), dimnames = list(rownames(mat), groups))
    for (i in seq_along(groups)) {
        group = groups[i] 
        col_group = mat[, group_vec == group, drop = FALSE]
        if (dim(col_group)[2] == 1) {
            final_mat[, i] = col_group
            which_mat[, i] = colnames(col_group)
        } else { 
            best_inds = apply(col_group, 1, which.max)
            best_queries = colnames(col_group)[best_inds]
            best_scores = col_group[cbind(1:num_conds, best_inds)]
            final_mat[, i] = best_scores
            which_mat[, i] = best_queries
        }
    }   
    
    return(list(best_scores = final_mat, best_queries = which_mat))
}

compute_per_gene_set_pvals_zscores = function(condition_x_gene_set_pred_sum_mat, control_condition_x_gene_set_pred_sum_mat_list, control_condition_x_gene_set_pred_sum_means_list, control_condition_x_gene_set_pred_sum_stdevs_list, control_types, alternative = c('greater', 'less', 'two-sided')) {
    # determine which direction of effect to detect
    alternative <- match.arg(alternative) 
    empirical_pval_fn <- switch(
        alternative,
        greater = empirical_pval_greater,
        less = empirical_pval_less,
        `two-sided` = empirical_pval_two_sided
    )
    
    # Get some dimensions
    num_control_types = length(control_types)
    num_conditions = dim(condition_x_gene_set_pred_sum_mat)[1]
    num_gene_sets = dim(condition_x_gene_set_pred_sum_mat)[2]

    ##### P-VALUES #####
    # Get p-values computed using different control types
    # as null distributions
    # Note: I compute the pvals for ALL conditions that I
    # have already specified I want to predict targets for,
    # not just for the treatments
    

    # Assemble output data structure!
    per_gene_set_pval_mats_by_control = replicate(num_control_types, 
                                                  matrix(numeric(num_conditions * num_gene_sets),
                                                         nrow = num_conditions,
                                                         dimnames = dimnames(condition_x_gene_set_pred_sum_mat)),
                                                  simplify = FALSE)
    names(per_gene_set_pval_mats_by_control) = control_types

    # Yes, this is a for loop in R. Satan's getting a bit chilled right now...
    for (i in 1:num_control_types) {

        for (j in 1:num_gene_sets) {
            
            cat(sprintf('computing %s-derived per-gene-set pval for gene set %s/%s\r', control_types[i], j, num_gene_sets))
            # message(sprintf('computing %s-derived per-gene-set pval for gene set %s/%s', control_types[i], j, num_gene_sets))
            per_gene_set_pval_mats_by_control[[i]][, j] = empirical_pval_fn(condition_x_gene_set_pred_sum_mat[, j],
                                                                            control_condition_x_gene_set_pred_sum_mat_list[[i]][, j])

        }
        cat('\n\n')

    }
    
    message('completed per-gene-set pval calculations')
    message(str(per_gene_set_pval_mats_by_control))

    ##### Z-SCORES #####
    # Get z-scores computed using different control types
    # as null distributions
    
    # Assemble output data structure!
    per_gene_set_zscore_mats_by_control = replicate(num_control_types, 
                                                    matrix(numeric(num_conditions * num_gene_sets),
                                                           nrow = num_conditions,
                                                           dimnames = dimnames(condition_x_gene_set_pred_sum_mat)),
                                                    simplify = FALSE)
    names(per_gene_set_zscore_mats_by_control) = control_types

    # Yes, this is a for loop in R. Satan's getting a bit chilled right now...
    for (i in 1:num_control_types) {

        for (j in 1:num_conditions) {
            
            cat(sprintf('computing %s-derived per-gene-set zscore for condition %s/%s\r', control_types[i], j, num_conditions))
            #message(sprintf('computing %s-derived per-gene-set zscore for condition %s/%s', control_types[i], j, num_conditions))
            per_gene_set_zscore_mats_by_control[[i]][j, ] = 
                    (condition_x_gene_set_pred_sum_mat[j, ] - control_condition_x_gene_set_pred_sum_means_list[[i]]) /
                     control_condition_x_gene_set_pred_sum_stdevs_list[[i]]
            
        }
        cat('\n\n')

    }
    
    message('Assembling final data structure...')

    # print(str(drug_subset_go_pred_sum_mat))

    per_gene_set_final_list = vector('list', num_control_types)
    names(per_gene_set_final_list) = control_types

    for (i in 1:num_control_types) {
        
        per_gene_set_final_list[[i]] = list(pval = per_gene_set_pval_mats_by_control[[i]],
                                            zscore = per_gene_set_zscore_mats_by_control[[i]])
    }

    return(per_gene_set_final_list)

}

compute_per_condition_pvals_zscores_2 = function(condition_x_gene_score_mat, condition_x_gene_set_sum_mat, gene_set_mat, num_per_cond_rand, gene_set_sizes, condition_x_gene_score_means, condition_x_gene_score_stdevs, seed, alternative = c('greater', 'less', 'two-sided')) {
    
    # determine which direction of effect to detect
    alternative = match.arg(alternative)
    
    # Get some dimensions
    num_conditions = dim(condition_x_gene_set_sum_mat)[1]
    num_gene_sets = dim(condition_x_gene_set_sum_mat)[2]
    num_genes = dim(condition_x_gene_score_mat)[2]
    
    ##### P-VALUES #####
    # Tricky part: compute p values for each per-drug
    # zscore.
    # Current scheme: shuffling gene labels on target_prediction_scores

    # Preallocate the output data structures before running
    # the loop that computes the # of randoms that beat real
    # predictions for each drug X gene set combination
    rand_result_mat = matrix(numeric(num_conditions * num_gene_sets),
                             nrow = num_conditions,
                             dimnames = dimnames(condition_x_gene_set_sum_mat))

    rand_beats_real_count_mat = matrix(integer(num_conditions * num_gene_sets),
                                       nrow = num_conditions,
                                       dimnames = dimnames(condition_x_gene_set_sum_mat))



    # Seed the random number generator right before
    # performing the randomizations!!!
    set.seed(seed)

    for (i in 1:num_per_cond_rand) {
        
        cat(sprintf('computing per-condition p values, randomization: %s/%s\r', i, num_per_cond_rand))
        #message(sprintf('computing per-condition p values, randomization: %s/%s', i, num_per_cond_rand))

        # Shuffle the gene set matrix (this is faster
        # than shuffling the condition X gene set
        # prediction matrix)
        shuffled_gene_set_mat = gene_set_mat[sample(x = 1:num_genes, size = num_genes, replace = FALSE), ]

        # Multiply the gene-level target prediction scores by
        # the shuffled gene set matrix to get per-condition
        # randomized sums! Assign the result to the same variable
        # that occupies the same position in memory to reduce
        # memory use issues (have not tested exactly how much
        # benefit comes from this - there is a slight speed
        # decrease)
        rand_result_mat[] = condition_x_gene_score_mat %*% shuffled_gene_set_mat

        # Determine which condition X gene set predictions using
        # the shuffled gene sets (shuffled on the gene side)
        # had the same or larger sum of gene-level target score
        # across all the genes in the gene set of interest.
        # Use subsetting to make sure the matrix is modified
        # in place to save memory
        rand_beats_real_count_mat[] = rand_beats_real_count_mat + as.integer(rand_result_mat >= condition_x_gene_set_sum_mat)

    }
    cat('\n\n')
    
    per_condition_pval_mat = switch(
        alternative,
        greater = rand_beats_real_count_mat / num_per_cond_rand,
        less = 1 - rand_beats_real_count_mat / num_per_cond_rand,
        `two-sided` = 2 * pmin(rand_beats_real_count_mat / num_per_cond_rand)
    )

    ##### Z-SCORES #####
    # Preallocate the output z-score matrix
    per_condition_zscore_mat = matrix(numeric(num_conditions * num_gene_sets),
                                      nrow = num_conditions,
                                      dimnames = dimnames(condition_x_gene_set_sum_mat))

    # Yes, another for loop. They work great lol.
    for (i in 1:num_gene_sets) {

        cat(sprintf('computing per-condition zscores for gene set %s/%s\r', i, num_gene_sets))
        #message(sprintf('computing per-condition zscores for gene set %s/%s', i, num_gene_sets))
        per_condition_zscore_mat[, i] = sqrt(gene_set_sizes[i]) * ((condition_x_gene_set_sum_mat[, i] / gene_set_sizes[i]) - condition_x_gene_score_means) / condition_x_gene_score_stdevs

    }
    cat('\n\n')

    per_condition_final_list <- list(pval = per_condition_pval_mat, zscore = per_condition_zscore_mat)
    
    return(per_condition_final_list)

}


# Define a function to extract the worst p value (or z score, if the
# p values match) between the two different control types (dmso and
# randomized drugs)
get_worst_case_pval_zscore <- function(scheme1_pval_mat, scheme1_zscore_mat, scheme2_pval_mat, scheme2_zscore_mat, scheme1_name, scheme2_name, alternative = c('greater', 'less', 'two-sided')) {

    alternative = match.arg(alternative)
    
    zscore_compare_fn <- switch(
        alternative,
        greater = `>`,
        less = `<`,
        `two-sided` = function(x,y) abs(x) > abs(y)
    )
    
    # Get sets of indices corresponding to where I should select
    # the scheme1 control-derived result vs the scheme2 drug-derived
    # result
    scheme1_beats_scheme2_pval <- scheme1_pval_mat < scheme2_pval_mat
    scheme2_beats_scheme1_pval <- scheme2_pval_mat < scheme1_pval_mat
    scheme1_equal_scheme2_pval <- scheme1_pval_mat == scheme2_pval_mat

    scheme1_beats_scheme2_zscore <- zscore_compare_fn(scheme1_zscore_mat, scheme2_zscore_mat)
    scheme2_beats_scheme1_zscore <- zscore_compare_fn(scheme2_zscore_mat, scheme1_zscore_mat)
    scheme1_equal_scheme2_zscore <- scheme1_zscore_mat == scheme2_zscore_mat

    pick_scheme1_indices <- scheme2_beats_scheme1_pval | (scheme1_equal_scheme2_pval & scheme2_beats_scheme1_zscore)
    pick_scheme2_indices <- scheme1_beats_scheme2_pval | (scheme1_equal_scheme2_pval & scheme1_beats_scheme2_zscore)
    tie_indices <- scheme1_equal_scheme2_pval & scheme1_equal_scheme2_zscore

    mat_dim <- dim(scheme1_pval_mat)
    mat_dimnames <- dimnames(scheme1_pval_mat)

    # Initialize some matrices
    worst_pval_mat <- matrix(numeric(), nrow = mat_dim[1], ncol = mat_dim[2], dimnames = mat_dimnames)
    worst_zscore_mat <- matrix(numeric(), nrow = mat_dim[1], ncol = mat_dim[2], dimnames = mat_dimnames)
    worst_control_name_mat <- matrix(character(), nrow = mat_dim[1], ncol = mat_dim[2], dimnames = mat_dimnames)

    # Fill in my worst case pval/zscore/control_name matrices!
    worst_pval_mat[pick_scheme1_indices] <- scheme1_pval_mat[pick_scheme1_indices]
    worst_pval_mat[pick_scheme2_indices] <- scheme2_pval_mat[pick_scheme2_indices]
    worst_pval_mat[tie_indices] <- scheme1_pval_mat[tie_indices]

    worst_zscore_mat[pick_scheme1_indices] <- scheme1_zscore_mat[pick_scheme1_indices]
    worst_zscore_mat[pick_scheme2_indices] <- scheme2_zscore_mat[pick_scheme2_indices]
    worst_zscore_mat[tie_indices] <- scheme1_zscore_mat[tie_indices]

    worst_control_name_mat[pick_scheme1_indices] <- scheme1_name
    worst_control_name_mat[pick_scheme2_indices] <- scheme2_name
    worst_control_name_mat[tie_indices] <- 'tie'

    list(worst_pval = worst_pval_mat,
         worst_zscore = worst_zscore_mat,
         control_name = worst_control_name_mat
         )
}

# Function to export diagnostic plots for any two schemes that generate a matrix of pvals
# and a matrix of zscores

print_diagnostic_plots_pval_zscore_2_schemes <- function(scheme1_pval_mat, scheme1_zscore_mat, scheme2_pval_mat, scheme2_zscore_mat, scheme1_name, scheme2_name, plots_folder) {

    dir.create(plots_folder, recursive = TRUE)

    master_plot_list <- list(
        # pval vs zscore plots
        scheme1_pval_vs_zscore_plot = qplot(x = as.vector(scheme1_zscore_mat),
                                             y = -log10(as.vector(scheme1_pval_mat)),
                                             xlab = sprintf('%s z-score', scheme1_name),
                                             ylab = sprintf('-log10 %s p value', scheme1_name),
                                             main = sprintf('Comparison of all %s-derived\nz-scores and p values', scheme1_name)
                                             ),

        scheme2_pval_vs_zscore_plot = qplot(x = as.vector(scheme2_zscore_mat),
                                             y = -log10(as.vector(scheme2_pval_mat)),
                                             xlab = sprintf('%s z-score', scheme2_name),
                                             ylab = sprintf('-log10 %s p value', scheme2_name),
                                             main = sprintf('Comparison of all %s-derived\nz-scores and p values', scheme2_name)
                                             ),

        # Across scheme pval vs pval and zscore vs zscore plots
        scheme1_vs_scheme2_pval_plot = qplot(x = -log10(as.vector(scheme1_pval_mat)),
                                              y = -log10(as.vector(scheme2_pval_mat)),
                                              xlab = sprintf('-log10 %s p value', scheme1_name),
                                              ylab = sprintf('-log10 %s p value', scheme2_name),
                                              main = sprintf('Comparison of all %s- and %s-derived p values', scheme1_name, scheme2_name)
                                              ),
        
        scheme1_vs_scheme2_zscore_plot = qplot(x = as.vector(scheme1_zscore_mat),
                                              y = as.vector(scheme2_zscore_mat),
                                              xlab = sprintf('%s z-score', scheme1_name),
                                              ylab = sprintf('%s z-score', scheme2_name),
                                              main = sprintf('Comparison of all %s- and %s-derived z-scores', scheme1_name, scheme2_name)
                                              ),

        # Plot only each condition's top p values or top zscores from each scheme
        # against each other
        scheme1_vs_scheme2_top_pval_plot = qplot(x = -log10(apply(scheme1_pval_mat, 1, min)),
                                              y = -log10(apply(scheme2_pval_mat, 1, min)),
                                              xlab = sprintf('-log10 %s min p value\nfor each condition', scheme1_name),
                                              ylab = sprintf('-log10 %s min p value\nfor each condition', scheme2_name),
                                              main = sprintf('Comparison of all %s- and %s-derived p values', scheme1_name, scheme2_name)
                                              ),

        scheme1_vs_scheme2_top_zscore_plot = qplot(x = apply(scheme1_zscore_mat, 1, max),
                                              y = apply(scheme2_zscore_mat, 1, max),
                                              xlab = sprintf('%s max z score\nfor each condition', scheme1_name),
                                              ylab = sprintf('%s max z score\nfor each condition', scheme2_name),
                                              main = sprintf('Comparison of all %s- and %s-derived z scores', scheme1_name, scheme2_name)
                                              )
    )

    # Add some common parameters to each plot for aesthetics

    # And some code to print all of these out to plots
    qc_plot <- function(ggplot_obj, filename) {
        png(file = file.path(plots_folder, filename), height = 4, width = 4, units = 'in', dpi = 200)
        print(ggplot_obj)
        dev.off()
    }

    # Change names of list to reflect actual identities
    new_names <- names(master_plot_list)
    new_names <- vapply(new_names, function(x) sub(pattern = 'scheme1', replacement = scheme1_name, x = x), character(1))
    new_names <- vapply(new_names, function(x) sub(pattern = 'scheme2', replacement = scheme2_name, x = x), character(1))
    names(master_plot_list) <- new_names

    lapply(names(master_plot_list), function(x) {
           qc_plot(master_plot_list[[x]], sprintf('%s.png', x))
    })
    
    return(NULL)

}


