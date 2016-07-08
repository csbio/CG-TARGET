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

get_gene_set_target_prediction_filename = function(folder, rand_scheme, num_rand, seed, similarity, per_cond_seed, per_cond_num_rand, gene_set_name, test_run) {
    file.path(folder, sprintf('gene_set_target_prediction_resampled-%s_%s-rands_seed-%s_%s_per-cond-seed-%s_%s-per-cond-rands_%s-gene-sets.txt.gz', rand_scheme, num_rand, seed, similarity,  per_cond_seed, per_cond_num_rand, gene_set_name, ifelse(test_run, 'test-run', '')))
}
