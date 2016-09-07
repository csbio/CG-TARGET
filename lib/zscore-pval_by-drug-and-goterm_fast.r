#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

library(foreach)
library(iterators)
library(ggplot2)
library(data.table)
library(Cairo)

# Source in libraries specific to this part of the script!
TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/one-sided_empirical_pval.R'))

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

compute_per_gene_set_pvals_zscores = function(condition_x_gene_set_pred_sum_mat, control_condition_x_gene_set_pred_sum_mat_list, control_condition_x_gene_set_pred_sum_means_list, control_condition_x_gene_set_pred_sum_stdevs_list, control_types) {

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
            per_gene_set_pval_mats_by_control[[i]][, j] = empirical_pval_greater(condition_x_gene_set_pred_sum_mat[, j],
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
            
            #message(sprintf('computing %s-derived per-gene-set pval for gene set %s/%s', control_types[i], j, num_gene_sets))
            #per_gene_set_pval_mats_by_control[[i]][, j] = empirical_pval_greater(condition_x_gene_set_pred_sum_mat[, j],
            #                                                                control_condition_x_gene_set_pred_sum_mat_list[[i]][, j])
        }
        cat('\n\n')

    }
    
    
    ## Get p values computed using different control types
    ## as null distributions
    ## Note: I compute the pvals for ALL conditions that I
    ## specify I want to predict targets for, not just for
    ## the treatments
    #per_go_pval_mats_by_control <- foreach(control_go_pred_sum_mat = control_go_pred_sum_mat_list, control_type = control_types) %do% {
    #    i_control_go_pred_sum_col <- iter(control_go_pred_sum_mat, by = 'col')
    #    i_drug_subset_go_pred_sum_col <- iter(drug_subset_go_pred_sum_mat, by = 'col')

    #    n = dim(control_go_pred_sum_mat)[2]

    #    # Note: previous bug here - I used dim(control_go_pred_sum_mat)[2], which
    #    # meant that all my pvalues were off by the ratio of (# of samples) to (# of
    #    # go terms)
    #    num_controls <- dim(control_go_pred_sum_mat)[1]
    #    per_go_pval_noexport = c('drug_gene_score_means', 'drug_gene_score_stdevs', 'go_term_sizes',
    #                             'drug_subset_go_pred_sum_mat', 'target_prediction_subset_mat',
    #                             'control_go_pred_sum_mat_list', 'control_go_pred_sum_stdevs_list',
    #                             'control_go_pred_sum_mat')
    #    per_go_pval_mat <- foreach(drug_subset_go_pred_sum_col = i_drug_subset_go_pred_sum_col, control_go_pred_sum_col = i_control_go_pred_sum_col, i = icount(), .combine = cbind, .maxcombine = 100000, .noexport = per_go_pval_noexport) %do% {

    #        message(sprintf('computing %s-derived per-gene-set pval for gene set %s/%s', control_type, i, n))

    #        drug_subset_go_pred_sum_col <- as.vector(drug_subset_go_pred_sum_col)
    #        control_go_pred_sum_col <- as.vector(control_go_pred_sum_col)

    #        # The other bug I had was that I counted the number of times
    #        # each treatment prediction beat the control predictions.
    #        # Instead, I want the number of times the control predictions
    #        # match or beat the control conditions.
    #        res = empirical_pval_greater(drug_subset_go_pred_sum_col, control_go_pred_sum_col)
    #        #res = vapply(drug_subset_go_pred_sum_col, function(x) {
    #        #       sum(x <= control_go_pred_sum_col) / num_controls
    #        #}, numeric(1))
    #        message(sprintf('%s/%s', i, n))
    #        message(str(res))
    #        #gc()
    #        res
    #    }
    #    message(str(drug_subset_go_pred_sum_mat))
    #    message(str(per_go_pval_mat))
    #    dimnames(per_go_pval_mat) <- dimnames(drug_subset_go_pred_sum_mat)
    #    per_go_pval_mat
    #}
    #
    #message('completed per-gene-set pval calculations')
    #message(str(per_gene_set_pval_mats_by_control))
    #
    # names(per_go_pval_mats_by_control) <- control_types

    # Get per-gene-set zscore calculations using different control
    # types as null distributions
    #per_go_zscore_mats_by_control <- foreach(control_go_pred_sum_means = control_go_pred_sum_means_list, control_go_pred_sum_stdevs = control_go_pred_sum_stdevs_list, control_type = control_types) %do% {
    #    i_drug_subset_go_pred_sum_row <- iter(drug_subset_go_pred_sum_mat, by = 'row')

    #    n = dim(drug_subset_go_pred_sum_mat)[1]
    #    per_go_zscore_mat <- foreach(drug_subset_go_pred_sum_row = i_drug_subset_go_pred_sum_row, i = icount(), .combine = rbind, .maxcombine = 100000) %dopar% {
    #        message(sprintf('computing %s-derived per-gene-set zscore for condition %s/%s', control_type, i, n))
    #        drug_subset_go_pred_sum_row <- as.vector(drug_subset_go_pred_sum_row)
    #        (drug_subset_go_pred_sum_row - control_go_pred_sum_means) / control_go_pred_sum_stdevs
    #    }
    #    dimnames(per_go_zscore_mat) <- dimnames(drug_subset_go_pred_sum_mat)
    #    per_go_zscore_mat
    #}

    #names(per_go_zscore_mats_by_control) <- control_types
    
    message('Assembling final data structure...')

    # print(str(drug_subset_go_pred_sum_mat))

    per_gene_set_final_list = vector('list', num_control_types)
    names(per_gene_set_final_list) = control_types

    for (i in 1:num_control_types) {
        
        per_gene_set_final_list[[i]] = list(pval = per_gene_set_pval_mats_by_control[[i]],
                                            zscore = per_gene_set_zscore_mats_by_control[[i]])
    }

    #per_go_final_list <- foreach(pval_mat = per_go_pval_mats_by_control, zscore_mat = per_go_zscore_mats_by_control) %do% {
    #    list(pval = pval_mat, zscore = zscore_mat)
    #}
    
    #names(per_go_final_list) <- control_types

    return(per_gene_set_final_list)

}

compute_per_condition_pvals_zscores = function(target_prediction_subset_mat, drug_subset_go_pred_sum_mat, num_by_drug_rand, go_term_sizes, drug_gene_score_means, drug_gene_score_stdevs, sample_type_split_vec, types_to_predict_for) {
    # Tricky part: compute p values for each per-drug
    # zscore.
    # Current scheme: shuffling gene labels on target_prediction_scores

    i_drug_subset_gene_score_row <- iter(target_prediction_subset_mat, by = 'row')
    i_drug_subset_go_pred_sum_row <- iter(drug_subset_go_pred_sum_mat, by = 'row')

    n = dim(drug_subset_go_pred_sum_mat)[1]

    per_drug_pval_mat <- foreach(drug_subset_gene_score_row = i_drug_subset_gene_score_row, drug_subset_go_pred_sum_row = i_drug_subset_go_pred_sum_row, i = icount(), .combine = rbind, .maxcombine = 100000) %dorng% {
        
        message(sprintf('computing per-condition p values for condition %s/%s', i, n))
        
        # Tidy up my data
        drug_subset_gene_score_row <- as.vector(drug_subset_gene_score_row)
        drug_subset_go_pred_sum_row <- as.vector(drug_subset_go_pred_sum_row)

        # Shuffle the gene labels on the drug target prediction vector
        num_genes <- length(drug_subset_gene_score_row)
        rand_drug_gene_score_mat <- do.call(rbind, replicate(num_by_drug_rand, sample(drug_subset_gene_score_row, size = num_genes, replace = FALSE), simplify = FALSE))

        # Obtain a distribution of target prediction score sums for each drug and GO term
        rand_go_pred_sum_mat <- rand_drug_gene_score_mat %*% go_term_mat

        # Compute p values by comparing the drug_subset_go_pred_sums to their corresponding random distributions
        # Recent edit: changed function to '>=', as p value definition is number of events
        # in the null distribution that are **as or more extreme** than the observation
        drug_subset_go_pred_sum_beats_rand_sum_mat <- sweep(rand_go_pred_sum_mat, MARGIN = 2, STATS = drug_subset_go_pred_sum_row, FUN = '>=')
        drug_subset_go_pred_pval <- colSums(drug_subset_go_pred_sum_beats_rand_sum_mat) / num_by_drug_rand
    }

    dimnames(per_drug_zscore_mat) <- dimnames(drug_subset_go_pred_sum_mat)
    dimnames(per_drug_pval_mat) <- dimnames(drug_subset_go_pred_sum_mat)
 


    # Compute per-condition z-scores
    i_drug_subset_go_pred_sum_col <- iter(drug_subset_go_pred_sum_mat, by = 'col')
    n = dim(drug_subset_go_pred_sum_mat)[2]
    per_drug_zscore_mat <- foreach(drug_subset_go_pred_sum_col = i_drug_subset_go_pred_sum_col, go_term_size = go_term_sizes, i = icount(), .combine = 'cbind', .maxcombine = 100000) %dopar% {
        message(sprintf('computing per-condition zscores for gene set %s/%s', i, n))
        drug_subset_go_pred_sum_col <- as.vector(drug_subset_go_pred_sum_col)

        drug_subset_gene_score_means <- drug_gene_score_means[sample_type_split_vec %in% types_to_predict_for]
        drug_subset_gene_score_stdevs <- drug_gene_score_stdevs[sample_type_split_vec %in% types_to_predict_for]

        sqrt(go_term_size) * ((drug_subset_go_pred_sum_col / go_term_size) - drug_subset_gene_score_means) / drug_subset_gene_score_stdevs
    }

    per_drug_final_list <- list(pval = per_drug_pval_mat, zscore = per_drug_zscore_mat)
    
    return(per_drug_final_list)

}

compute_per_condition_pvals_zscores_2 = function(condition_x_gene_score_mat, condition_x_gene_set_sum_mat, gene_set_mat, num_per_cond_rand, gene_set_sizes, condition_x_gene_score_means, condition_x_gene_score_stdevs, seed) {
    
    # Get some dimensions
    num_conditions = dim(condition_x_gene_set_sum_mat)[1]
    num_gene_sets = dim(condition_x_gene_set_sum_mat)[2]
    num_genes = dim(condition_x_gene_score_mat)[2]
    
    ##### P-VALUES #####
    # Tricky part: compute p values for each per-drug
    # zscore.
    # Current scheme: shuffling gene labels on target_prediction_scores

    #i_drug_subset_gene_score_row <- iter(target_prediction_subset_mat, by = 'row')
    #i_drug_subset_go_pred_sum_row <- iter(drug_subset_go_pred_sum_mat, by = 'row')

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
   
    per_condition_pval_mat = rand_beats_real_count_mat / num_per_cond_rand
    
    #per_drug_pval_mat <- foreach(drug_subset_gene_score_row = i_drug_subset_gene_score_row, drug_subset_go_pred_sum_row = i_drug_subset_go_pred_sum_row, i = icount(), .combine = rbind, .maxcombine = 100000) %dorng% {
    #    
    #    message(sprintf('computing per-condition p values for condition %s/%s', i, n))
    #    
    #    # Tidy up my data
    #    drug_subset_gene_score_row <- as.vector(drug_subset_gene_score_row)
    #    drug_subset_go_pred_sum_row <- as.vector(drug_subset_go_pred_sum_row)

    #    # Shuffle the gene labels on the drug target prediction vector
    #    num_genes <- length(drug_subset_gene_score_row)
    #    rand_drug_gene_score_mat <- do.call(rbind, replicate(num_by_drug_rand, sample(drug_subset_gene_score_row, size = num_genes, replace = FALSE), simplify = FALSE))

    #    # Obtain a distribution of target prediction score sums for each drug and GO term
    #    rand_go_pred_sum_mat <- rand_drug_gene_score_mat %*% go_term_mat

    #    # Compute p values by comparing the drug_subset_go_pred_sums to their corresponding random distributions
    #    # Recent edit: changed function to '>=', as p value definition is number of events
    #    # in the null distribution that are **as or more extreme** than the observation
    #    drug_subset_go_pred_sum_beats_rand_sum_mat <- sweep(rand_go_pred_sum_mat, MARGIN = 2, STATS = drug_subset_go_pred_sum_row, FUN = '>=')
    #    drug_subset_go_pred_pval <- colSums(drug_subset_go_pred_sum_beats_rand_sum_mat) / num_by_drug_rand
    #}
    #
    #dimnames(per_drug_pval_mat) <- dimnames(drug_subset_go_pred_sum_mat)
 

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


    ## Compute per-condition z-scores
    #i_condition_x_gene_set_sum_col <- iter(condition_x_gene_set_sum_mat, by = 'col')
    #n = dim(condition_x_gene_set_sum_mat)[2]
    #per_condition_zscore_mat <- foreach(condition_x_gene_set_sum_col = i_condition_x_gene_set_sum_col, gene_set_size = gene_set_sizes, i = icount(), .combine = 'cbind', .maxcombine = 100000) %dopar% {
    #    message(sprintf('computing per-condition zscores for gene set %s/%s', i, n))
    #    condition_x_gene_set_sum_col <- as.vector(condition_x_gene_set_sum_col)

    #    sqrt(gene_set_size) * ((condition_x_gene_set_sum_col / gene_set_size) - condition_x_gene_score_means) / condition_x_gene_score_stdevs
    #}
    #
    #dimnames(per_condition_zscore_mat) <- dimnames(condition_x_gene_set_sum_mat)

    per_condition_final_list <- list(pval = per_condition_pval_mat, zscore = per_condition_zscore_mat)
    
    return(per_condition_final_list)

}

compute_zscores_pvals_by_go_and_drug <- function(target_prediction_mat, go_term_mat, sample_type_split_vec, control_types, types_to_predict_for, num_by_drug_rand, load_point, save_points, gene_set_outdir) {

    # For the case of saving/loading partially processed data, define
    # an output folder!
    intermed_folder = file.path(gene_set_outdir, 'intermediate_data')
    dir.create(intermed_folder, recursive = TRUE)

    # Also, generate the potential filenames for saved intermediate
    # data
    save_1_filename = file.path(intermed_folder, 'save_point_1.RData')
    save_2_filename = file.path(intermed_folder, 'save_point_2.RData')

    # If the user does not specify a load point (or explicitly wants
    # to load from the beginning), then run the code from the
    # beginning!
    if (load_point < 1) {

        drug_gene_score_means <- rowMeans(target_prediction_mat)
        drug_gene_score_stdevs <- apply(target_prediction_mat, 1, sd)

        go_term_sizes <- colSums(go_term_mat)

        drug_go_pred_sum_mat <- target_prediction_mat %*% go_term_mat

        #     sample_types <- unique(sample_type_split_vec)

        #     control_types <- sample_types[sample_types != treatment_sample_type]
        
        drug_subset_go_pred_sum_mat <- drug_go_pred_sum_mat[sample_type_split_vec %in% types_to_predict_for, ]
        target_prediction_subset_mat <- target_prediction_mat[sample_type_split_vec %in% types_to_predict_for, ]
        
        # Change order of control_types so that the first control type has the
        # largest number of conditions. If they tie, order() will not change
        # their order.
        sample_type_counts = table(sample_type_split_vec)
        control_type_counts = sample_type_counts[control_types]
        control_order = order(control_type_counts, decreasing = TRUE)

        # Here's the final ordering step
        control_types = control_types[control_order]

        control_go_pred_sum_mat_list <- foreach(control_type = control_types) %do% {
            drug_go_pred_sum_mat[sample_type_split_vec %in% control_type, ]
        }


        rm(target_prediction_mat)
        rm(drug_go_pred_sum_mat)

        gc()

        control_go_pred_sum_means_list <- foreach(control_go_pred_sum_mat = control_go_pred_sum_mat_list) %do% {
            colMeans(control_go_pred_sum_mat)
        }
        
        control_go_pred_sum_stdevs_list <- foreach(control_go_pred_sum_mat = control_go_pred_sum_mat_list) %do% {
            apply(control_go_pred_sum_mat, 2, sd)
        }

        names(control_go_pred_sum_means_list) <- control_types
        names(control_go_pred_sum_stdevs_list) <- control_types

    }

    # End of code section 1

    # Here is very important testing/rerunning code. If the user wants to save
    # at "save point 1," then the R environment is saved here before moving on.
    # If the user wants to load the data previously saved at "save point 1,"
    # then the previously saved R environment is loaded here.
    if (1 %in% save_points) {
        save.image(save_1_filename)
    }

    if (load_point == 1) {
        load(save_1_filename)
    }

    # If the user-specified load point is less than 2, then this
    # next section of code must be run (because the data were
    # loaded in somewhere before this code and must be processed!)
    if (load_point < 2) {

        # Get p values computed using different control types
        # as null distributions
        # Note: I compute the pvals for ALL conditions, not just
        # the treatments
        per_go_pval_mats_by_control <- foreach(control_go_pred_sum_mat = control_go_pred_sum_mat_list, control_type = control_types) %do% {
            i_control_go_pred_sum_col <- iter(control_go_pred_sum_mat, by = 'col')
            i_drug_subset_go_pred_sum_col <- iter(drug_subset_go_pred_sum_mat, by = 'col')

            n = dim(control_go_pred_sum_mat)[2]
            
            # Note: previous bug here - I used dim(control_go_pred_sum_mat)[2], which
            # meant that all my pvalues were off by the ratio of (# of samples) to (# of
            # go terms)
            num_controls <- dim(control_go_pred_sum_mat)[1]
            per_go_pval_noexport = c('drug_gene_score_means', 'drug_gene_score_stdevs', 'go_term_sizes',
                                     'drug_subset_go_pred_sum_mat', 'target_prediction_subset_mat',
                                     'control_go_pred_sum_mat_list', 'control_go_pred_sum_stdevs_list',
                                     'control_go_pred_sum_mat')
            per_go_pval_mat <- foreach(drug_subset_go_pred_sum_col = i_drug_subset_go_pred_sum_col, control_go_pred_sum_col = i_control_go_pred_sum_col, i = icount(), .combine = cbind, .maxcombine = 100000, .noexport = per_go_pval_noexport) %dopar% {

                message(sprintf('computing %s-derived per-gene-set pval for gene set %s/%s', control_type, i, n))

                drug_subset_go_pred_sum_col <- as.vector(drug_subset_go_pred_sum_col)
                control_go_pred_sum_col <- as.vector(control_go_pred_sum_col)
               
                # The other bug I had was that I counted the number of times
                # each treatment prediction beat the control predictions.
                # Instead, I want the number of times the control predictions
                # match or beat the control conditions.
                res = vapply(drug_subset_go_pred_sum_col, function(x) {
                       sum(x <= control_go_pred_sum_col) / num_controls
                }, numeric(1))
                message(sprintf('%s/%s', i, n))
                message(str(res))
                gc()
                res
            }
            message(str(drug_subset_go_pred_sum_mat))
            message(str(per_go_pval_mat))
            dimnames(per_go_pval_mat) <- dimnames(drug_subset_go_pred_sum_mat)
            per_go_pval_mat
        }

        message('completed per-gene-set pval calculations')
        message(str(per_go_pval_mat))

        names(per_go_pval_mats_by_control) <- control_types

        # Get per-GO zscore calculations using different control
        # types as null distributions
        per_go_zscore_mats_by_control <- foreach(control_go_pred_sum_means = control_go_pred_sum_means_list, control_go_pred_sum_stdevs = control_go_pred_sum_stdevs_list, control_type = control_types) %do% {
            i_drug_subset_go_pred_sum_row <- iter(drug_subset_go_pred_sum_mat, by = 'row')

            n = dim(drug_subset_go_pred_sum_mat)[1]
            per_go_zscore_mat <- foreach(drug_subset_go_pred_sum_row = i_drug_subset_go_pred_sum_row, i = icount(), .combine = rbind, .maxcombine = 100000) %dopar% {
                message(sprintf('computing %s-derived per-gene-set zscore for condition %s/%s', control_type, i, n))
                drug_subset_go_pred_sum_row <- as.vector(drug_subset_go_pred_sum_row)
                (drug_subset_go_pred_sum_row - control_go_pred_sum_means) / control_go_pred_sum_stdevs
            }
            dimnames(per_go_zscore_mat) <- dimnames(drug_subset_go_pred_sum_mat)
            per_go_zscore_mat
        }

        names(per_go_zscore_mats_by_control) <- control_types

        # Compute per-drug z-scores
        i_drug_subset_go_pred_sum_col <- iter(drug_subset_go_pred_sum_mat, by = 'col')
        n = dim(drug_subset_go_pred_sum_mat)[2]
        per_drug_zscore_mat <- foreach(drug_subset_go_pred_sum_col = i_drug_subset_go_pred_sum_col, go_term_size = go_term_sizes, i = icount(), .combine = 'cbind', .maxcombine = 100000) %dopar% {
            message(sprintf('computing per-condition zscores for gene set %s/%s', i, n))
            drug_subset_go_pred_sum_col <- as.vector(drug_subset_go_pred_sum_col)

            drug_subset_gene_score_means <- drug_gene_score_means[sample_type_split_vec %in% types_to_predict_for]
            drug_subset_gene_score_stdevs <- drug_gene_score_stdevs[sample_type_split_vec %in% types_to_predict_for]

            sqrt(go_term_size) * ((drug_subset_go_pred_sum_col / go_term_size) - drug_subset_gene_score_means) / drug_subset_gene_score_stdevs
        }
    }

    # End of code section 2

    # More important testing/rerunning code. If the user wants to save
    # at "save point 2," then the R environment is saved here before moving on.
    # If the user wants to load the data previously saved at "save point 2,"
    # then the previously saved R environment is loaded here.
    if (2 %in% save_points) {
        save.image(save_2_filename)
    }

    if (load_point == 2) {
        load(save_2_filename)
    }

    # Run the following section of code if and only if the user is reading
    # the data in from files or loading in partially processed data before
    # this code is run.
    if (load_point < 3) {

        # Tricky part: compute p values for each per-drug
        # zscore.
        # Current scheme: shuffling gene labels on target_prediction_scores

        i_drug_subset_gene_score_row <- iter(target_prediction_subset_mat, by = 'row')
        i_drug_subset_go_pred_sum_row <- iter(drug_subset_go_pred_sum_mat, by = 'row')

        n = dim(drug_subset_go_pred_sum_mat)[1]

        per_drug_pval_mat <- foreach(drug_subset_gene_score_row = i_drug_subset_gene_score_row, drug_subset_go_pred_sum_row = i_drug_subset_go_pred_sum_row, i = icount(), .combine = rbind, .maxcombine = 100000) %dopar% {
            
            message(sprintf('computing per-condition p values for condition %s/%s', i, n))
            
            # Tidy up my data
            drug_subset_gene_score_row <- as.vector(drug_subset_gene_score_row)
            drug_subset_go_pred_sum_row <- as.vector(drug_subset_go_pred_sum_row)

            # Shuffle the gene labels on the drug target prediction vector
            num_genes <- length(drug_subset_gene_score_row)
            rand_drug_gene_score_mat <- do.call(rbind, replicate(num_by_drug_rand, sample(drug_subset_gene_score_row, size = num_genes, replace = FALSE), simplify = FALSE))

            # Obtain a distribution of target prediction score sums for each drug and GO term
            rand_go_pred_sum_mat <- rand_drug_gene_score_mat %*% go_term_mat

            # Compute p values by comparing the drug_subset_go_pred_sums to their corresponding random distributions
            # Recent edit: changed function to '>=', as p value definition is number of events
            # in the null distribution that are **as or more extreme** than the observation
            drug_subset_go_pred_sum_beats_rand_sum_mat <- sweep(rand_go_pred_sum_mat, MARGIN = 2, STATS = drug_subset_go_pred_sum_row, FUN = '>=')
            drug_subset_go_pred_pval <- colSums(drug_subset_go_pred_sum_beats_rand_sum_mat) / num_by_drug_rand
        }

        dimnames(per_drug_zscore_mat) <- dimnames(drug_subset_go_pred_sum_mat)
        dimnames(per_drug_pval_mat) <- dimnames(drug_subset_go_pred_sum_mat)

        
        # Assemble all of my target process predictions into a large
        # list data structure.

        message('Assembling final data structure...')

        print(str(drug_subset_go_pred_sum_mat))

        per_go_final_list <- foreach(pval_mat = per_go_pval_mats_by_control, zscore_mat = per_go_zscore_mats_by_control) %do% {
            list(pval = pval_mat, zscore = zscore_mat)
        }

        names(per_go_final_list) <- control_types

        per_drug_final_list <- list(pval = per_drug_pval_mat, zscore = per_drug_zscore_mat)

        list(per_gene_set = per_go_final_list, per_condition = per_drug_final_list)
    }

    # End of code section 3

}

# Define a function to extract the worst p value (or z score, if the
# p values match) between the two different control types (dmso and
# randomized drugs)
get_worst_case_pval_zscore <- function(scheme1_pval_mat, scheme1_zscore_mat, scheme2_pval_mat, scheme2_zscore_mat, scheme1_name, scheme2_name) {

    # Get sets of indices corresponding to where I should select
    # the scheme1 control-derived result vs the scheme2omized drug-derived
    # result
    scheme1_beats_scheme2_pval <- scheme1_pval_mat < scheme2_pval_mat
    scheme2_beats_scheme1_pval <- scheme1_pval_mat > scheme2_pval_mat
    scheme1_equal_scheme2_pval <- scheme1_pval_mat == scheme2_pval_mat

    scheme1_beats_scheme2_zscore <- scheme1_zscore_mat > scheme2_zscore_mat
    scheme2_beats_scheme1_zscore <- scheme1_zscore_mat < scheme2_zscore_mat
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
        CairoPNG(file = file.path(plots_folder, filename), height = 4, width = 4, units = 'in', dpi = 200)
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

get_sample_types <- function(drugnames, field) {

    vapply(strsplit(drugnames, '_'), `[`, character(1), field)
}



# Looks like I accidentally reimplemented this without knowing that I had written this...
# Oh well'-log10 %s p value', scheme1_name)
#
#get_process_prediction_drivers <- function(target_prediction_mat, go_term_mat, cutoff_expr) {
#    
#    cutoff_ex <- substitute(cutoff_expr)
#    
#    go_term_mat_log <- go_term_mat
#    mode(go_term_mat_log) <- 'logical'
#
#
#    i_go_mat_col <- iter(go_term_mat, by = 'col')
#    drivers_dts_list_by_go <- foreach(go_mat_col = i_go_mat_col, go = dimnames(go_term_mat)[[2]], .maxcombine = 100000) %dopar% {
#
#        tp_subset <- target_prediction_mat[, as.vector(go_mat_col)]
#
#        i_tp_subset_row <- iter(tp_subset, by = 'row')
#        drivers_dts_list <- foreach(x = i_tp_subset_row, drug = dimnames(tp_subset)[[1]], .maxcombine = 100000) %do% {
#
#            x <- as.vector(x)
#
#            ind <- eval(cutoff_expr)
#            
#            top_driver_vals <- x[ind]
#            top_driver_orfs <- names(tp_subset)[ind]
#
#            ordered_ind <- order(top_driver_vals, decreasing = TRUE)
#            ordered_top_driver_vals <- top_driver_vals[ordered_ind]
#            ordered_top_driver_orfs <- top_driver_orfs[ordered_ind]
#            ordered_top_driver_common <- orf_2_common(ordered_top_driver_orfs)
#
#            data.table(drug = drug,
#              GO = go,
#              score = paste(ordered_top_driver_vals, collapse = '|'),
#              orf = paste(ordered_top_driver_orfs, collapse = '|'),
#              common = paste(ordered_top_driver_common, collapse = '|')
#              )
#        }
#        rbindlist(drivers_dts_list)
#    }
#    rbindlist(drivers_dts_list_by_go)
#
#}
        






#
#    # Create plots folder and subdirectories, if they do not already exist
#    dir.create(plots_folder)
#    subdir_paths <- c(file.path(plots_folder, 'per-GO', c('within_controls', 'across_controls')),
#                      file.path(plots_folder, 'per-drug', 
#
#
#    # Get the top prediction for each drug, determined first by smallest
#    # p value and then largest z-score (if the p values are both 0,
#    # for example)
#    top_drug_prediction_by_control <- foreach(pval_mat = per_go_pval_mats_by_control, zscore_mat = per_go_zscore_mats_by_control) %do% {
#        drug_pvals <- iter(pval_mat, by = 'row')
#        drug_zscores <- iter(zscore_mats, by = 'row')
#        foreach(drug_pvals = i_drug_pvals, drug_zscore = i_drug_zscores, combine = 'rbind') %do% {
#            top_index <- order(drug_pvals, -drug_zscore)[1]
#            c(drug_pvals[top_index], drug_zscore[top_index])
#        }
#    }
#
#
#    # Export plots of each type of p value against each other,
#    # all types of z-score against each other, and p value vs.
#    # zscore within each control type.
#    for (i in length(control_types)) {
#
#        # Plot -log p value vs zscore for each control type
#
#
#
#        for (j in (i + 1):length(control_types) {
#
#             # Plot -log10 of pvalues against each other
#             neg_log10_pvals_vec_x <- -log10(as.vector(per_go_pval_mats_by_control[[i]]))
#             neg_log10_pvals_vec_y <- -log10(as.vector(per_go_pval_mats_by_control[[j]]))
#             pval_plot <- qplot(x = neg_log10_pvals_vec_x,
#                                y = neg_log10_pvals_vec_y,
#                                main = sprintf('%s-derived vs. %s-derived\np values', control_types[j], control_types[i])
#                                xlab = sprintf('-log10(%s-derived p value)', control_types[i])
#                                ylab = sprintf('-log10(%s-derived p value)', control_types[j])
#                                )
#             CairoPNG(
#
#             # Plot z-scores against each other
#
#
#
#    # Get a per-go p value matrix representing the maximum (worst)
#    # p value between dmso-based p value and the randomized sample-
#    # based p value.
#

