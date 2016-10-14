#!/usr/bin/env Rscript

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

## CHANGE THIS HEADING DOCUMENTATION!!!!!!!!
#
#
## This script will take in the results of chemical genomics target prediction
## and return "z scores" describing the relationship between a drug's average
## target prediction score for genes that match a GO term, and the drug's overall
## average target prediction score.
#
## This script exports all zscore calculations
## Analysis is done using other scripts
#


library(optparse)
library(yaml)
library(reshape2)
library(data.table)

# Source in libraries for gene set target prediction
TARGET_PATH = Sys.getenv('TARGET_PATH')
source(file.path(TARGET_PATH, 'lib/go_drivers.r'))
source(file.path(TARGET_PATH, 'lib/file_naming.r'))
source(file.path(TARGET_PATH, 'lib/top_table.r'))
source(file.path(TARGET_PATH, 'lib/table_to_globular.r'))
source(file.path(TARGET_PATH, 'lib/zscore-pval_by-drug-and-goterm_fast.r'))
source(file.path(TARGET_PATH, 'lib/pos_args.r'))
source(file.path(TARGET_PATH, 'lib/filenames.r'))
source(file.path(TARGET_PATH, 'lib/datasets.r'))


# print(sessionInfo())

# Some function definitions...
save_data_factory = function(folder) {
    # This function returns a function that saves data in a
    # folder that is specified at the time of function creation.
    # The resulting function only requires a "save_point"
    # argument, which is x.
    
    force(folder)
    fun = function(x) {
        outfile = file.path(folder, sprintf('save_point_%s.RData', x))
        message(sprintf('\nsaving data at save point %s...', x))
        message(sprintf('location: %s\n', outfile))
        # Get a list of everything except for functions, and also
        # remove the save_points and load_point variablex. We want
        # both functions and save/load points to be defined for
        # each run of the code. (don't want old bugs hanging
        # around or wasting time on saving when the user doesn't
        # want it)
        # If I think about it, removing the load_point variable
        # probably doesn't matter, but I'm preventing it from
        # saving anyway.
        object_vec = setdiff(ls(envir = .GlobalEnv), lsf.str(envir = .GlobalEnv))
        object_vec = setdiff(object_vec, c('save_points', 'load_point'))
        #print(object_vec)
        save(list = object_vec, file = outfile, envir = .GlobalEnv)
        message('done')
    }
    return(fun)
}


load_data_factory = function(folder) {
    # This function returns a function that saves loads data
    # from a folder that is specified at the time of function
    # creation. The resulting function only requires a
    # "load_point" argument, which is x.
    
    force(folder)
    fun = function(x) {
        infile = file.path(folder, sprintf('save_point_%s.RData', x))
        message(sprintf('\nloading data from save point %s...', x))
        message(sprintf('location: %s\n', infile))
        load(infile, envir = globalenv())
        message('done')
    }
    return(fun)
}


# Parse the arguments, as these need to be read regardless of which starting
# point in the code the user wants to start the analysis from.
# Note: while it is totally ok (and more or less the point) for the "load_point"
# and "save_points" arguments to change between runs on the same dataset
# (for example, stepping through the analysis, one save point at a time),
# "--test_run" only matters when no load point is specified (or it is 0).
positional_arguments_list = c(CONFIG_FILE = 'yaml-formatted gene-set prediction configuration file.')

option_list = list(make_option('--test_run', action = 'store_true', default = FALSE, help = 'Performs predictions for a small subset of the data, to ensure that everything is working properly before the full set of predictions are made (can take days for large screens). This option only has an effect when "--load_point" is not specified or is 0. If "--test_run" has been specified and the analysis saved halfway through using "--save_points," then the analysis will still be in "test run" mode when loading that save point, and vice versa if "--test_run" was not specified. NOTE: the results from this run will not be accurate!!!'),
                   make_option('--load_point', default = 0, type = 'integer', help = 'Point in the processing at which a previously partially-processed, saved dataset should be loaded and its processing continued through the rest of the script. Possible values are in c(0, 1, 2). 0 (default): run all code from the beginning; 1: resume the analysis at per-gene-set p-value computations; 2: resume the analysis at per-condition p-value computations. While primarily for testing purposes, this could be useful during real analyses.'),
                   make_option('--save_points', default = '', type = 'character', help = 'Point(s) in the processing at which to save the data in their present, partially-processed form. See "--load_point" for descriptions of these points. If the user wishes to specify multiple save points, these must be comma-delimited *without* spaces (e.g. --save_points 1,2). Data are saved in "*.RData" format.')
                   )
final_usage_positional = get_usage_and_positional(positional_arguments_list)

parser = OptionParser(usage = final_usage_positional, option_list = option_list)
arguments = parse_args(parser, positional_arguments = length(positional_arguments_list))
opt = arguments$options
arg = arguments$args

print(opt)
print(arg)



####### First, deal with the save/load point arguments
# Load 'em in...
load_point = opt$load_point
save_points_raw = opt$save_points
save_points_split = strsplit(save_points_raw, ',')[[1]]
save_points = as.integer(save_points_split)

# Make sure no saving is supposed to happen before loading
if (any(save_points < load_point)) {
    stop('Specified intermediate file save points must occur after the specified load point!')
}

# Looks like I have to parse the config file, *just a bit*, so
# that it knows where the output directory is. Grrrrrrr...
config_f = file(arg[1], 'rt')
config_params = yaml.load_file(config_f)
close(config_f)

# Get the output folder for intermediate data
# (and create if it doesn't exist!)
output_table_folder = get_gene_set_target_prediction_folder(config_params$Required_arguments$output_folder)
intermed_folder = file.path(output_table_folder, 'intermediate_data')
dir.create(intermed_folder, recursive = TRUE)

# Generate the save_data and load_data functions
save_data = save_data_factory(intermed_folder)
load_data = load_data_factory(intermed_folder)

# If the user does not specify a load point (or explicitly wants
# to load from the beginning), then run the code from the
# beginning! Otherwise, all of this code will be skipped, as the
# options have already been defined and saved in the user's
# R environment from the first round of processing!

# Now, on to the real stuff!
if (load_point < 1) {


    #positional_arguments_list = c(gene_predictions = 'Gzipped file containing the gene-level target predictions for each condition.',
    #                              rand_gene_predictions = 'Gzipped file containing the gene-level target predictions for each resampled condition.',
    #                              gene_sets = 'Two-column, tab-delimited file (WITH HEADER) that maps gene identifiers (must have same column name as corresponding column in the query info table) to gene sets.',
    #                              query_info_table = 'Table with a "query_key" column that maps to the query identifiers in the genetic interaction data and another column that maps to the gene identifiers in the gene_sets file.',
    #                              driver_cutoff = 'Similarity value (similarity to the condition profile) above which a gene is reported as a "driver" of particular gene-set prediction. Default value cannot be specified because unbounded similarity measures are accepted as input.',
    #                              output_table = 'File into which the results are written. (it is gzipped, so the extension should be *.txt.gz)'
    #                              )

    #option_list = list(
    #    make_option('--sample_table', help = 'File with a table that maps conditions to control information'),
    #    make_option('--neg_control_column', help = 'Name of True/False column in the sample_table that indicates which samples are to be used as negative controls in gene-set prediction. If this option is not specified, only the resampled conditions will be used for gene-set prediction.'),
    #    make_option('--condition_name_column', help = 'Name of a column in the sample table with human-readable condition names'),
    #    make_option(c('-n', '--per_condition_randomizations'), type = 'integer', help = 'The number of "per-condition" randomizations to perform. When not specified, the default is to perform the same number of randomizations that are present in the rand_gene_predictions file.'),
    #    make_option('--gene_set_names', help = 'Two-column, tab-delimited file that maps each gene set identifier to an interpretable name for gene-set prediciton output. The values in the first column must match the gene sets in the gene_sets file above. The values in the second column will be included in the final results table (along with those in the first column).'),
    #    make_option('--gene_name_column', help = 'Column in "query_info_table" that contains interpretable gene names. In S. cerevisiae, for example this would be the column that contains "TUB3" in the row for ORF YML124C'),
    #    make_option('--test_run', action = 'store_true', default = FALSE, help = 'Performs predictions for a small subset of the data, to ensure that everything is working properly before the full set of predictions are made (can take days for large screens). The results from this run will not be accurate!!!'),
    #    make_option('--num_cores', type = 'integer', default = 1, help = 'The number of cores used to run the gene set predictions!'),
    #    make_option(c('-v', '--verbosity'), type = 'integer', default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1, is default.')
    #    )


    ################# Load in the required files, deal with arguments

    # Note: I had concerns that not parsing the config file within
    # this block of code would lead the user to be able to alter
    # the config params on a saved run that was loaded halfway into
    # the analysis. I now believe this is unfounded, as all variables
    # will be restored to those from the config file at the time the
    # analysis was initiated when the user loads all variables in
    # from a specified save point.

    ####### Now for the parameters in the config file!!!
    output_folder = config_params$Required_arguments$output_folder
    gene_pred_folder = get_gene_target_folder(output_folder)
    gene_pred_file = get_gene_target_prediction_filename(gene_pred_folder,
                                                         config_params$Required_arguments$cg_gi_similarity_measure
                                                         )
    gene_prediction_tab = fread(sprintf('gzip -dc %s', gene_pred_file), header = TRUE, colClasses = c('character', 'character', 'character', 'numeric'))

    # Since the user can specify not to use the resampled profiles derived
    # from the entire dataset, here I start dealing with that. The code
    # will look a bit clunky, but this will be more reliable/predictable
    # than the alternative of creating an empty data.table to funnel
    # through all of the same steps. That will likely wreak havoc.
    zero_resampled_profiles = (config_params$Required_arguments$`per-array_resampling_scheme` == 0) | (config_params$Required_arguments$`num_per-array_resampled_profiles`== 0)

    if (!zero_resampled_profiles) {

        rand_gene_pred_file = get_gene_target_prediction_resampled_filename(
                                  gene_pred_folder,
                                  config_params$Required_arguments$`per-array_resampling_scheme`,
                                  config_params$Required_arguments$`num_per-array_resampled_profiles`,
                                  config_params$Required_arguments$`per-array_resampling_seed`,
                                  config_params$Required_arguments$cg_gi_similarity_measure
                                  )
        rand_gene_prediction_tab = fread(sprintf('gzip -dc %s', rand_gene_pred_file), header = TRUE, colClasses = c('character', 'character', 'character', 'numeric'))
    } else {

        # If there are no resampled profiles, make sure that the config_params,
        # which are used later in the script, accurately reflect that.
        config_params$Required_arguments$`per-array_resampling_scheme` = 0
        config_params$Required_arguments$`num_per-array_resampled_profiles` = 0
        config_params$Required_arguments$`per-array_resampling_seed` = 0
    }

    sample_table_filename = config_params$Required_arguments$cg_col_info_table
    sample_table = fread(sample_table_filename, header = TRUE, colClasses = 'character')

    gene_set_name = config_params$Required_arguments$gene_set_name
    gi_data_name = config_params$Required_arguments$gi_dataset_name
    gene_set_info = get_gene_set_info(gene_set_name, TARGET_PATH)
    print(gene_set_info)
    gene_set_tab_full = fread(gene_set_info$filename, header = TRUE, sep = '\t', colClasses = 'character')
    gene_set_gene_id_col = names(gene_set_tab_full)[1]
    gene_set_id_col = names(gene_set_tab_full)[2]
    gene_set_name_col = gene_set_info$interpretable_column
    print(gene_set_tab_full)
    print(gene_set_id_col)
    print(gene_set_name_col)
    gene_set_tab = gene_set_tab_full[, c(gene_set_gene_id_col, gene_set_id_col), with = FALSE]

    gi_info = get_gi_info(config_params$Required_arguments$gi_dataset_name, TARGET_PATH)
    #gi_tab = fread(sprintf('gzip -dc %s', gi_filenames$gi_tab))
    query_info_tab = fread(gi_info$gi_query_tab)

    #query_info_file = config_params$Required_arguments$gi_query_info_table
    #query_info_tab = fread(query_info_file, header = TRUE, colClasses = 'character')

    driver_cutoff = as.numeric(config_params$Required_arguments$driver_cutoff)
    if (is.na(driver_cutoff)) { stop('driver_cutoff argument must be numeric.') }

    output_table_folder = get_gene_set_target_prediction_folder(output_folder)
    output_table = get_gene_set_target_prediction_filename(
                       output_table_folder,
                       config_params$Required_arguments$`per-array_resampling_scheme`,
                       config_params$Required_arguments$`num_per-array_resampled_profiles`,
                       config_params$Required_arguments$`per-array_resampling_seed`,
                       config_params$Required_arguments$cg_gi_similarity_measure,
                       config_params$Required_arguments$`per-condition_randomization_seed`,
                       config_params$Required_arguments$`num_per-condition_randomizations`,
                       config_params$Required_arguments$gene_set_name,
                       opt$test_run
                       )

    dir.create(output_table_folder, recursive = TRUE)

    true_false_control_map = c('TRUE' = 'expt_control', 'FALSE' = 'treatment', 'True' = 'expt_control', 'False' = 'treatment')

    ################### Load in optional files and options
    # Gene set ID to interpretable name table
    # Might want to change this to "if (!gene_set_name_col == '')" instead of expecting it to be null or NA
    if (!is.null(gene_set_name_col)) {
        #gene_set_name_file = config_params$Options$gene_set_target_prediction$gene_set_name_table
        gene_set_name_tab = unique(gene_set_tab_full[, c(gene_set_id_col, gene_set_name_col), with = FALSE])
        print(gene_set_name_tab)
        setnames(gene_set_name_tab, c('gene_set', 'gene_set_name'))
    } else {
        gene_set_name_tab = NULL
    }

    # Note: if no column is specified containing interpretable gene names,
    # then this variable will already be NULL.
    gene_name_column = gi_info$query_genename_col
    if (!is.null(gene_name_column)) {
        if (! gene_name_column %in% names(query_info_tab)) {
            stop(sprintf('column %s was specified as the column in the query info table that contains interpretable gene names, but the columns does not exist! query info table is here:\n%s', gene_name_column, gi_info$gi_query_tab))
        }
    }

    ################ Read in the sample table and negative control/condition name columns, if they exist
    # Handle the negative control column
    neg_control_col = config_params$Options$gene_set_target_prediction$negative_control_column
    if (!is.null(neg_control_col)) {
        if (! neg_control_col %in% names(sample_table)) {
            stop(sprintf('column %s not in cg_col_info_table!', neg_control_col))
        }
        control_map = true_false_control_map[sample_table[[neg_control_col]]]
        names(control_map) = sample_table[, sprintf('%s_%s', screen_name, expt_id)]
    } else {
        control_map = rep('treatment', dim(sample_table)[1])
        names(control_map) = sample_table[, sprintf('%s_%s', screen_name, expt_id)]
    }

    # Handle the condition name column
    cond_name_col = config_params$Options$gene_set_target_prediction$condition_name_column
    if (!is.null(cond_name_col)) {
        ### May need to revisit if the columns are split or not...
        if (! cond_name_col %in% names(sample_table)) {
            stop(sprintf('column %s not in cg_col_info_table!', cond_name_col))
        }
        condition_name_tab = sample_table[, list(condition = sprintf('%s_%s', screen_name, expt_id), condition_name = sample_table[[cond_name_col]])]
        
        # The condition_tab needs to have rand-by-strain conditions added so the join
        # doesn't cause us to lose all or the rand-by-strain conditions!
        # **But only if there is a table of rand-by-strain (resampled) gene target
        # predictions!**
        if (!zero_resampled_profiles) {
            rand_condition_name_tab = unique(rand_gene_prediction_tab[, list(screen_name, expt_id)])[, list(condition = sprintf('%s_%s', screen_name, expt_id), condition_name = sprintf('%s_%s', screen_name, expt_id))]
            condition_name_tab = rbind(condition_name_tab, rand_condition_name_tab)
        }

        print(condition_name_tab)
    } else {
        condition_name_tab = NULL
    }

    
    print(sample_table)
    print(control_map)
   
    ############ Here, convert all the tables into matrices and them remove
    ############ before any workspaces can be saved (lots of unnecessary
    ############ space would be taken up!

    ###### Smash some data around
    gene_prediction_mat = acast(gene_prediction_tab, screen_name + expt_id ~ query_key, value.var = 'score')
    if (!zero_resampled_profiles) {
        rand_gene_prediction_mat = acast(rand_gene_prediction_tab, screen_name + expt_id ~ query_key, value.var = 'score')
        all_prediction_mat = rbind(gene_prediction_mat, rand_gene_prediction_mat)
    
        ###### Line up the control vectors to the matrix with all predictions
        rand_control_map = rep('rand-by-strain', dim(unique(rand_gene_prediction_tab[, list(screen_name, expt_id)]))[1])
        names(rand_control_map) = unique(rand_gene_prediction_tab[, sprintf('%s_%s', screen_name, expt_id)])
        
        # Don't need these!
        rm(rand_gene_prediction_tab)
        rm(rand_gene_prediction_mat)

    } else {
        all_prediction_mat = gene_prediction_mat
    }

    # Don't need these!
    rm(gene_prediction_tab)
    rm(gene_prediction_mat)

    #all_prediction_tab = rbind(gene_prediction_tab, rand_gene_prediction_tab)
    #all_prediction_mat = acast(all_prediction_tab, screen_name + expt_id ~ query_key, value.var = 'score')


    gc()

}

# Here is very important testing/rerunning code. If the user wants to save
# at "save point 1," then the R environment is saved here before moving on.
# If the user wants to load the data previously saved at "save point 1,"
# then the previously saved R environment is loaded here.
if (1 %in% save_points) {
    save_data(1)
}
if (load_point == 1) {
    load_data(1)
}

# If the user-specified load point is less than 2, then this
# next section of code must be run (because the data were
# loaded in somewhere before this code and must be processed
# to completion!)
if (load_point < 2) {


    # Check for real and random conditions overlapping (should NEVER happen)
    # Only check if there are "random" (really, "resampled") conditions
    if (!zero_resampled_profiles) {
        print(control_map[1:10])
        print(rand_control_map[1:10])
        overlap_conds = intersect(names(control_map), names(rand_control_map))
        if (length(overlap_conds) > 0) {
            stop(sprintf('The following real and random conditions have the same name:\n%s', paste(overlap_conds, collapse = '\n')))
        }
    }

    # Combine the control maps and get one master control type vector!!!
    if (!zero_resampled_profiles) {
    
        all_control_map = c(control_map, rand_control_map)
    } else {
        all_control_map = control_map
    }

    all_control_vec = all_control_map[rownames(all_prediction_mat)]
    print(all_control_vec)
    print(sprintf('Number of %s conditions: %s', names(table(all_control_vec)), table(all_control_vec)))

    # If any sample types are still NA, filter them out of the dataset
    inds_to_remove = is.na(all_control_vec)
    conds_to_remove = rownames(all_prediction_mat)[inds_to_remove]

    all_control_vec = all_control_vec[!is.na(all_control_vec)]
    all_prediction_mat = all_prediction_mat[!is.na(all_control_vec), ]

    # Get a matrix with one column per unique query gene (there may be multiple alleles for the same gene),
    # and for each row (condition) in that column, we take the maximum prediction score against all
    # alleles of the same query gene. Also, get a matrix that indicates which allele the selected
    # max score came from.

    # First, get the vector of query genes that define where the multiple-allele cases are.
    #print(gene_set_tab)
    #print(query_info_tab)
    # WHY DO I DRAW THE NAME OF THE QUERY COLUMN FROM THE GENE_SET_TAB?!
    # Well now I don't :)
    query_column = gi_info$query_sys_name_col
    query_map = query_info_tab[[query_column]]
    names(query_map) = query_info_tab[['query_key']]
    query_gene_vec = query_map[colnames(all_prediction_mat)]

    print(str(all_prediction_mat))
    print(str(query_gene_vec))
    print(all_prediction_mat[1:5, (dim(all_prediction_mat)[2]-9):(dim(all_prediction_mat)[2])])
    print(query_gene_vec[(length(query_gene_vec) - 9):length(query_gene_vec)])

    # Get the matrix with the best predictions among multiple alleles
    best_query_mat_list = best_score_per_col_group(all_prediction_mat, query_gene_vec)

    best_prediction_mat = best_query_mat_list[['best_scores']]
    best_query_mat = best_query_mat_list[['best_queries']]

    # Get a map from query key column to an interpretable name, if the
    # gene_name_column exists. Then create a matrix that matches the
    # best_query_mat matrix and contains interpretable query gene
    # names. If that column was not given, create a matrix of NA
    # values instead.
    if (!is.null(gene_name_column)) {
        gene_name_map = query_info_tab[[gene_name_column]]
        names(gene_name_map) = query_info_tab[['query_key']]
        best_query_name_mat = matrix(gene_name_map[best_query_mat], nrow = nrow(best_query_mat), ncol = ncol(best_query_mat), dimnames = dimnames(best_query_mat))
    } else {
        best_query_name_mat = matrix(NA_character_, nrow = nrow(best_query_mat), ncol = ncol(best_query_mat), dimnames = dimnames(best_query_mat))
    }

    ###### Determine which sample types exist and will be used to make predictions
    # At some point, should the user be able to specify this instead of just assuming
    # that everything needs to have its target predicted?
    sample_types_to_predict_for <- unique(all_control_vec)

    #########################
    #########################

    #########
    ### What are the control types? (There may be experimental controls, or maybe not!)
    controls <- intersect(sample_types_to_predict_for, c('expt_control', 'rand-by-strain'))

    # Reshape the gene set table into a matrix
    gene_set_mat_formula = as.formula(sprintf('%s ~ %s', names(gene_set_tab)[1], names(gene_set_tab)[2]))
    gene_set_tab[, annotation := 1]
    gene_set_matrix = acast(data = gene_set_tab, formula = gene_set_mat_formula, value.var = 'annotation', fill = 0)

    ## Get a matrix of all possible propagated orf to GO term mappings
    #raw_go_matrix <- get_gene_go_bp_is_a_propagated_matrix(raw_orfs)
    #
    ## Only select GO terms with 4 to 200 genes from the lightnet
    ## mapped to them
    #go_term_sizes <- apply(raw_go_matrix, 2, sum)
    #go_matrix <- raw_go_matrix[, go_term_sizes >= 4 & go_term_sizes <= 200]

    # Match up the gene set matrix with the prediction/query matrices:
    common_genes = intersect(rownames(gene_set_matrix), colnames(best_prediction_mat))
    gene_set_matrix = gene_set_matrix[common_genes, , drop = FALSE]
    best_prediction_mat <- best_prediction_mat[, common_genes, drop = FALSE]
    best_query_mat <- best_query_mat[, common_genes, drop = FALSE]
    best_query_name_mat = best_query_name_mat[, common_genes, drop = FALSE]

    # Remove GO terms from the gene set matrix without any annotations (must only be
    # annotated from genes not in the set that is predicted against in the gene-level
    # target prediction step).
    nonzero_degree_gene_set_inds = colSums(gene_set_matrix) > 0
    gene_set_matrix = gene_set_matrix[, nonzero_degree_gene_set_inds, drop = FALSE]
    print(sprintf('Removed %s gene sets with no annotations', sum(!nonzero_degree_gene_set_inds)))


    # Get p values and zscores for each drug --> go combination
    # This will be done on a per-GO term and a per-drug basis, and
    # the per-GO computations can be performed using as many
    # different control types of samples as possible.  From this
    # function call, one gets back all p values and z scores for
    # ALL different computation methods.

    # First, create a vector that splits the matrix up into
    # treatments and controls
    # treat_control_vec <- vector('character', dim(tp_mat)[1])
    # treat_control_vec[] <- 'Treatment'
    # treat_control_vec[grepl('DMSO', dimnames(tp_mat)[[1]], fixed = TRUE)] <- 'DMSO'
    # treat_control_vec[grepl('Rand-by-strain', dimnames(tp_mat)[[1]], fixed = TRUE)] <- 'Rand-by-strain'


    # test_set <- c(which(treat_control_vec == 'Treatment')[1:20],
    #               which(treat_control_vec == 'DMSO')[1:20],
    #               which(treat_control_vec == 'Rand-by-strain')[1:20]
    #               )

    # Select a subset for debugging purposes
    # Selects 100 of each sample type (or fewer, if fewer than
    # 100 of a sample type
    test_inds <- do.call(c, lapply(unique(all_control_vec), function(x) {
                                   inds <- which(all_control_vec == x)
                                   inds[1:(min(length(inds), 100))]
    }))


    # Is this a test to make sure everything works, or is this the real deal?
    # Look for the command line option "--test_run"
    if (opt$test_run) {
        condition_inds <- test_inds
    } else {
        condition_inds <- TRUE
    }


    # stop('time to figure out what is going terribly wrong')

    ###### The following steps of code, until the end of this section, are
    ###### preparation for performing the z-score and p-value computations
    target_prediction_mat = best_prediction_mat[condition_inds,]
    sample_type_split_vec = all_control_vec[condition_inds]

    condition_gene_score_means <- rowMeans(target_prediction_mat)
    condition_gene_score_stdevs <- apply(target_prediction_mat, 1, sd)

    gene_set_sizes <- colSums(gene_set_matrix)

    condition_gene_set_pred_sum_mat <- target_prediction_mat %*% gene_set_matrix

    #     sample_types <- unique(sample_type_split_vec)

    #     controls <- sample_types[sample_types != treatment_sample_type]

    # Perform some subsetting to only include sample types to make
    # predictions for
    condition_subset_gene_set_pred_sum_mat <- condition_gene_set_pred_sum_mat[sample_type_split_vec %in% sample_types_to_predict_for, ]
    target_prediction_subset_mat <- target_prediction_mat[sample_type_split_vec %in% sample_types_to_predict_for, ]
    condition_subset_gene_score_means = condition_gene_score_means[sample_type_split_vec %in% sample_types_to_predict_for]
    condition_subset_gene_score_stdevs = condition_gene_score_stdevs[sample_type_split_vec %in% sample_types_to_predict_for]

    # Change order of controls so that the first control type has the
    # largest number of conditions. If they tie, order() will not change
    # their order.
    sample_type_counts = table(sample_type_split_vec)
    control_counts = sample_type_counts[controls]
    control_order = order(control_counts, decreasing = TRUE)

    # Here's the final ordering step
    controls = controls[control_order]

    #control_gene_set_pred_sum_mat_list <- foreach(control = controls) %do% {
    #    condition_gene_set_pred_sum_mat[sample_type_split_vec %in% control, ]
    #}

    control_gene_set_pred_sum_mat_list = vector('list', length(controls))
    for (i in seq_along(controls)) {
        control = controls[[i]]
        control_gene_set_pred_sum_mat_list[[i]] = condition_gene_set_pred_sum_mat[sample_type_split_vec == control, ]
    }
    names(control_gene_set_pred_sum_mat_list) = controls


    rm(target_prediction_mat)
    rm(condition_gene_set_pred_sum_mat)

    gc()

    #control_gene_set_pred_sum_means_list <- foreach(control_gene_set_pred_sum_mat = control_gene_set_pred_sum_mat_list) %do% {
    #    colMeans(control_gene_set_pred_sum_mat)
    #}

    #control_gene_set_pred_sum_stdevs_list <- foreach(control_gene_set_pred_sum_mat = control_gene_set_pred_sum_mat_list) %do% {
    #    apply(control_gene_set_pred_sum_mat, 2, sd)
    #}

    control_gene_set_pred_sum_means_list = lapply(control_gene_set_pred_sum_mat_list, colMeans)
    control_gene_set_pred_sum_stdevs_list = lapply(control_gene_set_pred_sum_mat_list,
                                                   function(x) apply(x, 2, sd))

    names(control_gene_set_pred_sum_means_list) = controls
    names(control_gene_set_pred_sum_stdevs_list) = controls

}


# Here is very important testing/rerunning code. If the user wants to save
# at "save point 2," then the R environment is saved here before moving on.
# If the user wants to load the data previously saved at "save point 2,"
# then the previously saved R environment is loaded here.
if (2 %in% save_points) {
    save_data(2)
}
if (load_point == 2) {
    load_data(2)
}


########### Get p-values and z-scores for per-gene-set evaluations

# If the user-specified load point is less than 3, then this
# next section of code must be run (because the data were
# loaded in somewhere before this code and must be processed
# to completion!)
if (load_point < 3) {

    if (!(is.null(neg_control_col) & zero_resampled_profiles)) {

        # Only compute per-gene-set scores if negative control and/or
        # resampled profiles are present. Otherwise, skip and replace
        # with a reasonable data structure.
        per_gene_set_pvals_zscores = compute_per_gene_set_pvals_zscores(condition_subset_gene_set_pred_sum_mat, control_gene_set_pred_sum_mat_list, control_gene_set_pred_sum_means_list, control_gene_set_pred_sum_stdevs_list, controls)
    } else {

        # Else, this list is just NULL instead. To be dealt with later.
        per_gene_set_pvals_zscores = NULL
    }

}

# Again, here is important testing/rerunning code. If the user wants to save
# at "save point 3," then the R environment is saved here before moving on.
# If the user wants to load the data previously saved at "save point 3,"
# then the previously saved R environment is loaded here.
if (3 %in% save_points) {
    save_data(3)
}
if (load_point == 3) {
    load_data(3)
}

########### Get p-values and z-scores for per-condition evaluations

# If the user-specified load point is less than 4, then this
# next section of code must be run (because the data were
# loaded in somewhere before this code and must be processed
# to completion!)
if (load_point < 4) {
    
    # Settle the per-condition randomization question
    num_per_condition_rand = config_params$Required_arguments$`num_per-condition_randomizations`

    # Get the seed (required input from the user)
    # Seed isn't set until right before randomizations
    # are performed, inside the function below.
    seed = config_params$Required_arguments$`per-condition_randomization_seed`
    print(seed)
    if (seed == 'rand') {
        seed = NULL
    } else if (!is.numeric(seed)) {
        stop('specified per-condition_randomization_seed is neither numeric nor "rand"')
    }
        
    per_condition_pvals_zscores = compute_per_condition_pvals_zscores_2(target_prediction_subset_mat, condition_subset_gene_set_pred_sum_mat, gene_set_matrix, num_per_condition_rand, gene_set_sizes, condition_subset_gene_score_means, condition_subset_gene_score_stdevs, seed)

    # Combine results into a final list and remove the individual lists
    all_pvals_zscores = list(per_gene_set = per_gene_set_pvals_zscores, per_condition = per_condition_pvals_zscores)

    # Remove any components of the list that are NULL. If both
    # are NULL, stop the code with an error since we need at
    # least one set of p-values and z-scores to move on.
    all_pvals_zscores = all_pvals_zscores[!vapply(all_pvals_zscores, is.null, logical(1))]

    if (length(all_pvals_zscores) == 0) {

        stop('There are no computed p-values/z-scores, from either the per-gene-set or per-condition analyses, with which to move forward. Please check your config file and your data to ensure you either have generated resampled profiles and/or have enough negative control profiles for per-gene-set analyses, and/or you have set the number of per-condition randomizations to a value greater than zero.')
    }

    if(!is.null(per_gene_set_pvals_zscores)) rm(per_gene_set_pvals_zscores)
    if(!is.null(per_condition_pvals_zscores)) rm(per_condition_pvals_zscores)
    gc()

}

# Save/load point 4!
if (4 %in% save_points) {
    save_data(4)
}
if (load_point == 4) {
    load_data(4)
}

# Now, get all versions of p values and zscores!
#system.time(all_pvals_zscores <- compute_zscores_pvals_by_go_and_drug(target_prediction_mat = best_prediction_mat[condition_inds,],
#                                                                      go_term_mat = gene_set_matrix,
#                                                                      sample_type_split_vec = all_control_vec[condition_inds],
#                                                                      control_types = controls,
#                                                                      types_to_predict_for = sample_types_to_predict_for,
#                                                                      num_per_condition_rand = num_per_condition_rand,
#                                                                      load_point = load_point,
#                                                                      save_points = save_points,
#                                                                      gene_set_outdir = output_table_folder))


######### Combine the p-values and z-scores so that there is one p-value
######### and z-score per condition X gene set pair. The largest p-value and its
######### accompanying z-score are chosen, to minimize the chance that the
######### results we call significant are actually easily observable among any
######### appropriate null models.

# If the data are already loaded, proceed with processing!
if (load_point < 5) {

    print(str(all_pvals_zscores))
    print(controls)

    #print(which(is.na(all_pvals_zscores[['per_gene_set']][[controls[1]]][['pval']]), arr.ind = TRUE))
    #print(sum(is.na(all_pvals_zscores[['per_gene_set']][[controls[1]]][['pval']])))

    #print(which(is.na(all_pvals_zscores[['per_gene_set']][[controls[1]]][['zscore']]), arr.ind = TRUE))
    #print(sum(is.na(all_pvals_zscores[['per_gene_set']][[controls[1]]][['zscore']])))

    #bad_go_term_inds = unique(which(is.na(all_pvals_zscores[['per_gene_set']][[controls[1]]][['zscore']]), arr.ind = TRUE)[, 2])
    #print(str(gene_set_matrix[, bad_go_term_inds]))
    #print(colSums(gene_set_matrix[, bad_go_term_inds]))

    #print(which(is.na(all_pvals_zscores[['per_gene_set']][[controls[2]]][['pval']])))
    #print(which(is.na(all_pvals_zscores[['per_gene_set']][[controls[2]]][['zscore']])))
    # get worst per-go pvals/zscores
    if (length(all_pvals_zscores[['per_gene_set']]) == 2) {
        per_gene_set_worst_case_mats <- get_worst_case_pval_zscore(all_pvals_zscores[['per_gene_set']][[controls[1]]][['pval']],
                                                                   all_pvals_zscores[['per_gene_set']][[controls[1]]][['zscore']],
                                                                   all_pvals_zscores[['per_gene_set']][[controls[2]]][['pval']],
                                                                   all_pvals_zscores[['per_gene_set']][[controls[2]]][['zscore']],
                                                                   controls[1],
                                                                   controls[2]
                                                                   )
    } else if (length(all_pvals_zscores[['per_gene_set']]) == 1) {
        print(str(all_pvals_zscores))
        per_gene_set_worst_case_mats = list(worst_pval = all_pvals_zscores[['per_gene_set']][[controls[1]]][['pval']],
                                            worst_zscore = all_pvals_zscores[['per_gene_set']][[controls[1]]][['zscore']],
                                            control_name = matrix(controls[1], nrow = dim(all_pvals_zscores[['per_gene_set']][[controls[1]]][['pval']])[1], ncol = dim(all_pvals_zscores[['per_gene_set']][[controls[1]]][['pval']])[2], dimnames = dimnames(all_pvals_zscores[['per_gene_set']][[controls[1]]]))
                                            )
    } else if (length(all_pvals_zscores[['per_gene_set']]) > 2) {

        # Since this code can now tolerate having no per_gene_set-
        # derived predictions, it's only a problem if there are
        # too many control types - more than there should be.
        stop(sprintf('Not sure how this happened, as there are %s control types when there should be either 1 or 2', length(controls)))
    }

    # get worst overall pvals/zscores (worst between
    # per-gene-set and per-condition, but ONLY if
    # there are both per_gene_set pvals/zscores and
    # per_condition pvals/zscores)
    if (length(all_pvals_zscores[['per_gene_set']]) == 0) {

        # If there are no per-gene-set predictions, then we just use
        # the per-condition predictions.
        overall_worst_case_mats = list(worst_pval = all_pvals_zscores[['per_condition']][['pval']],
                                       worst_zscore = all_pvals_zscores[['per_condition']][['zscore']],
                                       control_name = matrix('per_condition', nrow = dim(all_pvals_zscores[['per_condition']][['pval']])[1], ncol = dim(all_pvals_zscores[['per_condition']][['pval']])[2], dimnames = dimnames(all_pvals_zscores[['per_condition']][['pval']]))
                                       )

    } else if (length(all_pvals_zscores[['per_condition']]) == 0) {
        
        # Same thing for per-condition predictions
        overall_worst_case_mats = list(worst_pval = per_gene_set_worst_case_mats[['worst_pval']],
                                       worst_zscore = per_gene_set_worst_case_mats[['worst_zscore']],
                                       control_name = matrix('per_gene_set', nrow = dim(all_pvals_zscores[['per_gene_set']][['pval']])[1], ncol = dim(all_pvals_zscores[['per_gene_set']][['pval']])[2], dimnames = dimnames(all_pvals_zscores[['per_gene_set']][['pval']]))
                                       )
    } else {

        # Otherwise, there is at least one per-condition and one per-gene-set
        # set of statistics, this is how to combine them.
        overall_worst_case_mats <- get_worst_case_pval_zscore(per_gene_set_worst_case_mats$worst_pval,
                                                              per_gene_set_worst_case_mats$worst_zscore,
                                                              all_pvals_zscores$per_condition$pval,
                                                              all_pvals_zscores$per_condition$zscore,
                                                              'per_gene_set',
                                                              'per_condition'
                                                              )
    }

    # Specify which per-gene-set scheme was used if the worst pval/zscore
    # was generated using the per-gene-set scheme. But if there are no
    # per_gene_set indices, then don't do it!
    if (exists('per_gene_set_worst_case_mats')) {
        per_gene_set_indices <- overall_worst_case_mats$control_name == 'per_gene_set'
        overall_worst_case_mats$control_name[per_gene_set_indices] <- per_gene_set_worst_case_mats$control_name[per_gene_set_indices]
    }
}

# Save/load point 5!
if (5 %in% save_points) {
    save_data(5)
}
if (load_point == 5) {
    load_data(5)
}

########### Compute and assemble the drivers of the gene-set predictions

# If the data are already loaded, proceed with processing!
if (load_point < 6) {

    # I need to subset the target prediction matrix so I can get the correct
    # process prediction drivers back!
    best_prediction_subset_mat = best_prediction_mat[condition_inds, ][all_control_vec[condition_inds] %in% sample_types_to_predict_for, ]
    best_query_subset_mat = best_query_mat[condition_inds, ][all_control_vec[condition_inds] %in% sample_types_to_predict_for, ]
    best_query_name_subset_mat = best_query_name_mat[condition_inds, ][all_control_vec[condition_inds] %in% sample_types_to_predict_for, ]

    # Here I will get a table of drivers of my go process predictions
    gene_set_drivers_dt = get_gene_set_drivers(best_prediction_subset_mat, best_query_subset_mat, best_query_name_subset_mat, gene_set_matrix, cutoff = driver_cutoff)
    # If the driver gene names were not specified, remove that column from the target prediction table!
    if (is.null(gene_name_column)) {
        gene_set_drivers_dt[, driver_name := NULL]
    }

}

# Save/load point 6!
if (6 %in% save_points) {
    save_data(6)
}
if (load_point == 6) {
    load_data(6)
}

## Shut down MPI cluster, if it exists
#if (ncores > 1) {
#    closeCluster(clus)
#}

# print(gene_set_drivers_dt[condition == 'SHANGHAI-1511_000088'])

########### Table mashing and exporting!

# If the data are already loaded, proceed with processing!
if (load_point < 7) {

    # Set the key on the go_drivers_dt so that it can be joined with the
    # process prediction data.tables
    setkeyv(gene_set_drivers_dt, c('condition', 'gene_set'))

    ## Melt per-gene-set only predictions into a data table
    #per_gene_set_worst_case_df <- lapply(per_gene_set_worst_case_mats, as.vector)
    #per_gene_set_worst_case_df <- as.data.frame(per_gene_set_worst_case_df)
    #
    #per_gene_set_GENE_SET <- as.character(col(per_gene_set_worst_case_mats$worst_pval, as.factor = TRUE))
    #per_gene_set_CONDITION <- as.character(row(per_gene_set_worst_case_mats$worst_pval, as.factor = TRUE))
    #
    #per_go_worst_case_dt <- data.table(condition = per_gene_set_CONDITION, gene_set = per_gene_set_GENE_SET, per_gene_set_worst_case_df)

    # Melt overall only predictions into a data table
    overall_worst_case_df <- lapply(overall_worst_case_mats, as.vector)
    overall_worst_case_df <- as.data.frame(overall_worst_case_df)

    overall_GENE_SET <- as.character(col(overall_worst_case_mats$worst_pval, as.factor = TRUE))
    overall_CONDITION <- as.character(row(overall_worst_case_mats$worst_pval, as.factor = TRUE))

    overall_worst_case_dt <- data.table(condition = overall_CONDITION, gene_set = overall_GENE_SET, overall_worst_case_df)


    # Join the process score tables with the process driver table
    #setkeyv(per_go_worst_case_dt, c('condition', 'gene_set'))
    setkeyv(overall_worst_case_dt, c('condition', 'gene_set'))

    #per_gene_set_worst_case_dt <- gene_set_drivers_dt[per_gene_set_worst_case_dt]
    overall_worst_case_dt <- gene_set_drivers_dt[overall_worst_case_dt]

    # Add a column indicating control/treatment category for each condition.
    overall_worst_case_dt[, sample_type := all_control_vec[condition]]

    # If the user provided a table mapping gene set IDs to interpretable gene set names,
    # then add those to the final table!
    if (!is.null(gene_set_name_tab)) {
    #    setkey(per_gene_set_worst_case_dt, gene_set)
        setkey(overall_worst_case_dt, gene_set)
        setkey(gene_set_name_tab, gene_set)
    #    per_gene_set_worst_case_dt = per_gene_set_worst_case_dt[gene_set_name_tab]
        overall_worst_case_dt = overall_worst_case_dt[gene_set_name_tab, nomatch = 0]
    }

    # If the user provided a table mapping condition IDs to to interpretable condition names,
    # then add those to the final table!
    if (!is.null(condition_name_tab)) {
    #    setkey(per_gene_set_worst_case_dt, gene_set)
        setkey(overall_worst_case_dt, condition)
        setkey(condition_name_tab, condition)
    #    per_gene_set_worst_case_dt = per_gene_set_worst_case_dt[gene_set_name_tab]
        overall_worst_case_dt = overall_worst_case_dt[condition_name_tab, nomatch = 0]
    }


    ## Check size distribution of these go terms
    #go_terms <- per_go_worst_case_dt[, unique(GO)]
    #go_termsizes <- colSums(go_matrix[, go_terms])
    #hist(go_termsizes, breaks = 200)
    ## Comment: I think the size distribution checks out just fine :)
    ## However, these terms may be slightly different from the ones
    ## I have predicted in the past.  I am more confident in these
    ## propagations, however

    ##### Print out some diagnostic plots describing how all of the pvals/zscores
    ##### relate to each other across the types of control conditions I used to
    ##### derive them
    ####
    ##### Need to fix speed, or DO NOT RUN UNLESS YOU HAVE EXTRA TIME!!!
    ##### Also, plots probably look fugly at this point
    #####with(all_pvals_zscores$per_go, print_diagnostic_plots_pval_zscore_2_schemes(DMSO$pval,
    #####                                                                            DMSO$zscore,
    #####                                                                            `Rand-by-strain`$pval,
    #####                                                                            `Rand-by-strain`$zscore,
    #####                                                                            'DMSO',
    #####                                                                            'Rand-by-strain',
    #####                                                                            file.path(prefix, sprintf('%s_process-prediction_qc-plots', Sys.Date()))
    #####                                                                            ))
    #####
    #####print_diagnostic_plots_pval_zscore_2_schemes(per_go_worst_case_mats$worst_pval,
    #####                                             per_go_worst_case_mats$worst_zscore,
    #####                                             overall_worst_case_mats$worst_pval,
    #####                                             overall_worst_case_mats$worst_zscore,
    #####                                             'per-GO',
    #####                                             'per-drug',
    #####                                             file.path(prefix, sprintf('%s_process-prediction_qc-plots', Sys.Date()))
    #####                                             )
    #####
    #####
    #
    #
    ## Sort each prediction table by condition, pval, and zscore
    #setkeyv(per_go_worst_case_dt, c('drug'))
    #per_go_worst_case_dt <- per_go_worst_case_dt[, .SD[order(-worst_pval, worst_zscore, decreasing = TRUE)], by = list(drug)]
    #
    setkey(overall_worst_case_dt, condition)
    # overall_worst_case_dt <- overall_worst_case_dt[, .SD[order(-worst_pval, worst_zscore, decreasing = TRUE)], by = list(drug)]
    overall_worst_case_dt <- overall_worst_case_dt[order(worst_pval, -worst_zscore)]

    # Set the order of the columns!!! (This varies depending on which other tables were provided and joined to the predictions)
    print(names(overall_worst_case_dt))
    cols = c('condition', 'condition_name', 'gene_set', 'gene_set_name', 'p_value', 'z_score', 'driver_query', 'driver_name', 'driver_score', 'sample_type', 'control_name')
    setnames(overall_worst_case_dt, c('worst_pval', 'worst_zscore'), c('p_value', 'z_score'))
    if (is.null(gene_set_name_tab)) cols = c(cols[1:(which(cols == 'gene_set_name') - 1)], cols[(which(cols == 'gene_set_name') + 1):length(cols)])
    if (is.null(condition_name_tab)) cols = c(cols[1:(which(cols == 'condition_name') - 1)], cols[(which(cols == 'condition_name') + 1):length(cols)])
    if (is.null(gene_name_column)) cols = c(cols[1:(which(cols == 'driver_name') - 1)], cols[(which(cols == 'driver_name') + 1):length(cols)])
    setcolorder(overall_worst_case_dt, cols)

    # Write out to (compressed) table-formatted file
    #outfile_per_go <- gzfile(file.path(prefix, sprintf('%s_per-go_pval-zscore.txt.gz', Sys.Date())), open = 'wt')
    #write.table(per_go_worst_case_dt, file = outfile_per_go, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
    #close(outfile_per_go)


    outfile_overall <- gzfile(output_table, open = 'wb')
    write.table(overall_worst_case_dt, file = outfile_overall, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
    close(outfile_overall)

    # Write out a top table and a globular table for human consumption

    # First, make filenames for the top table and glob top table outfiles
    #out_top_per_go <- file.path(prefix, sprintf('%s_per-go_pval-zscore_top-ten.txt.gz', Sys.Date()))
    #out_top_glob_per_go <- file.path(prefix, sprintf('%s_per-go_pval-zscore_top-ten_glob.txt.gz', Sys.Date()))
    ## Now export prediction tables for just per_go predictions
    #top_ten_tab_per_go <- export_top_table(dat = per_go_worst_case_dt, outfile = out_top_per_go, select_expr = rank(worst_pval, ties.method = 'min') <= 10, split_by_expr = list(drug), order_within_by_expr = order(-worst_pval, worst_zscore, decreasing = TRUE), connection_FUN = 'gzfile')
    #table_to_glob(dat = top_ten_tab_per_go, outfile = out_top_glob_per_go, by_vec = c('drug'), connection_FUN = 'gzfile')

    # More filenames
    #out_top_overall <- file.path(prefix, sprintf('%s_overall_pval-zscore_top-ten.txt.gz', Sys.Date()))
    #out_top_glob_overall <- file.path(prefix, sprintf('%s_overall_pval-zscore_top-ten_glob.txt.gz', Sys.Date()))
    ## And, export the overall pval/zscore top table and glob top table files
    #top_ten_tab_overall <- export_top_table(dat = overall_worst_case_dt, outfile = out_top_overall, select_expr = rank(worst_pval, ties.method = 'min') <= 10, split_by_expr = list(drug), order_within_by_expr = order(-worst_pval, worst_zscore, decreasing = TRUE), connection_FUN = 'gzfile')
    #table_to_glob(dat = top_ten_tab_overall, outfile = out_top_glob_overall, by_vec = c('drug'), connection_FUN = 'gzfile')

}

# Save/load point 7! This would be purely for debugging, I think.
if (7 %in% save_points) {
    save_data(7)
}
if (load_point == 7) {
    load_data(7)
}



