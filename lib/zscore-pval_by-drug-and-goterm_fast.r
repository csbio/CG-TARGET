library(foreach)
library(iterators)
library(ggplot2)
library(data.table)
library(Cairo)

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


compute_zscores_pvals_by_go_and_drug <- function(target_prediction_mat, go_term_mat, sample_type_split_vec, control_types, types_to_predict_for, num_by_drug_rand) {

    drug_gene_score_means <- rowMeans(target_prediction_mat)
    drug_gene_score_stdevs <- apply(target_prediction_mat, 1, sd)

    go_term_sizes <- colSums(go_term_mat)

    drug_go_pred_sum_mat <- target_prediction_mat %*% go_term_mat

    #     sample_types <- unique(sample_type_split_vec)

    #     control_types <- sample_types[sample_types != treatment_sample_type]
    
    drug_subset_go_pred_sum_mat <- drug_go_pred_sum_mat[sample_type_split_vec %in% types_to_predict_for, ]
    target_prediction_subset_mat <- target_prediction_mat[sample_type_split_vec %in% types_to_predict_for, ]

    control_go_pred_sum_mat_list <- foreach(control_type = control_types) %do% {
        drug_go_pred_sum_mat[sample_type_split_vec %in% control_type, ]
    }

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
        per_go_pval_mat <- foreach(drug_subset_go_pred_sum_col = i_drug_subset_go_pred_sum_col, control_go_pred_sum_col = i_control_go_pred_sum_col, i = icount(), .combine = cbind, .maxcombine = 100000) %dopar% {

            message(sprintf('computing %s-derived per-gene-set pval for gene set %s/%s', control_type, i, n))

            drug_subset_go_pred_sum_col <- as.vector(drug_subset_go_pred_sum_col)
            control_go_pred_sum_col <- as.vector(control_go_pred_sum_col)
           
            # The other bug I had was that I counted the number of times
            # each treatment prediction beat the control predictions.
            # Instead, I want the number of times the control predictions
            # match or beat the control conditions.
            vapply(drug_subset_go_pred_sum_col, function(x) {
                   sum(x <= control_go_pred_sum_col) / num_controls
            }, numeric(1))
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

