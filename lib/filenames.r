#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# The functions in this R file are to unify file naming across
# all target prediction scripts.
get_resampled_profile_folder = function(outdir) {
    file.path(outdir, 'resampled_profiles')
}

get_resampled_profile_filename = function(folder, rand_scheme, num_rand, seed) {
    file.path(folder, sprintf('resampled_profiles_scheme-%s_%s-rands_seed-%s.txt.gz', rand_scheme, num_rand, seed))
}

get_gene_target_folder = function(outdir) {
    file.path(outdir, 'gene_target_prediction')
}

get_gene_target_prediction_filename = function(folder, similarity) {
    file.path(folder, sprintf('gene_target_prediction_%s.txt.gz', similarity))
}

get_gene_target_prediction_resampled_filename = function(folder, rand_scheme, num_rand, seed, similarity) {
    file.path(folder, sprintf('gene_target_prediction_resampled-%s_%s-rands_seed-%s_%s.txt.gz', rand_scheme, num_rand, seed, similarity))
}

get_gene_set_target_prediction_folder = function(outdir) {
    file.path(outdir, 'gene_set_target_prediction')
}

get_gene_set_target_prediction_filename_v1 = function(folder, rand_scheme, num_rand, seed, similarity, per_cond_seed, per_cond_num_rand, gene_set_name, test_run) {
    # If either the randomization scheme or the number of randomizations
    # is set to zero, then set both to zero, including the seed!
    if (rand_scheme == 0 | num_rand == 0) {
        rand_scheme = 0
        num_rand = 0
        seed = 0
    }
    file.path(folder, sprintf('gene_set_target_prediction_resampled-%s_%s-rands_seed-%s_%s_per-cond-seed-%s_%s-per-cond-rands_%s-gene-sets.txt.gz', rand_scheme, num_rand, seed, similarity,  per_cond_seed, per_cond_num_rand, gene_set_name, ifelse(test_run, 'test-run', '')))
}

get_gene_set_target_prediction_filename_v2 = function(folder, rand_scheme, num_rand, seed, similarity, per_cond_seed, per_cond_num_rand, gene_set_name, min_termsize, max_termsize, test_run) {
    # If either the randomization scheme or the number of randomizations
    # is set to zero, then set both to zero, including the seed!
    if (rand_scheme == 0 | num_rand == 0) {
        rand_scheme = 0
        num_rand = 0
        seed = 0
    }
    file.path(folder, sprintf('gene_set_target_prediction_%s_%s_%s_%s_%s_%s_%s_%s-%s%s.txt.gz', rand_scheme, num_rand, seed, similarity, per_cond_seed, per_cond_num_rand, gene_set_name, min_termsize, max_termsize, ifelse(test_run, '_test-run', '')))
}

get_final_results_folder = function(outdir, gene_set_name) {
    file.path(outdir, 'final_results', gene_set_name)
}

get_final_results_tabs_folder = function(outdir, gene_set_name) {
    final_results_folder = get_final_results_folder(outdir, gene_set_name)
    file.path(final_results_folder, 'tables')
}

get_final_results_plots_folder = function(outdir, gene_set_name) {
    final_results_folder = get_final_results_folder(outdir, gene_set_name)
    file.path(final_results_folder, 'plots')
}
