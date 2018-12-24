#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

library(data.table)
library(ggplot2)
library(reshape2)

source(file.path(TARGET_PATH, 'lib/scott_themes.r'))
source(file.path(TARGET_PATH, 'lib/table_printing.r'))

get_top_prediction_per_condition<- function(dat) {

    setkey(dat, condition)

    dat[, .SD[order(p_value, -z_score)][1], by = condition]
}

add_top_discovery_counts<- function(top_dat) {

    setkeyv(top_dat, 'sample_type')

    # Determine which control types were used
    control_types = intersect(c('expt_control', 'rand-by-strain'), unique(top_dat[['sample_type']]))
    
    # Within each sample type sort by ascending pval, 
    # then descending zscore. Then, add cumsum column
    top_dat <- top_dat[, .SD[order(p_value, -z_score)], by = sample_type]
    top_dat[, condition_count := as.numeric(1:.N), by = sample_type]

    top_dat[, scaled_condition_count := condition_count]
    top_dat[, normalization_factor := 1]
    for (x in control_types) {
        norm_factor_colname = sprintf('%s_normalization_factor', x)
        top_dat[[norm_factor_colname]] = top_dat[, sum(sample_type == 'treatment') / sum(sample_type == x)]
        top_dat[x, normalization_factor := top_dat[x][[norm_factor_colname]]]
    }
    top_dat[, scaled_condition_count := condition_count * normalization_factor]
    for (x in control_types) {
        norm_factor_colname = sprintf('%s_normalization_factor', x)
        top_dat[[norm_factor_colname]] = NULL
    }

#    top_dat[, dmso_normalization_factor := sum(sample_type == 'Treatment') / sum(sample_type == 'DMSO')]
#    top_dat[, rand_normalization_factor := sum(sample_type == 'Treatment') / sum(sample_type == 'Rand-by-strain')]
#    top_dat[, normalization_factor := 1]
#    top_dat['DMSO', normalization_factor := dmso_normalization_factor]
#    top_dat['Rand-by-strain', normalization_factor := rand_normalization_factor]
#    top_dat[, scaled_condition_count := condition_count * normalization_factor]
#    top_dat[, c('dmso_normalization_factor', 'rand_normalization_factor') := NULL]

    gc()

    top_dat
}

plot_discoveries_vs_pval <- function(top_counts_dt, filename_meat, pre) {

    # Do the plot thing    
    dvp_plot <- ggplot(data = top_counts_dt) +
        geom_line(aes(x = p_value, y = scaled_condition_count, color = sample_type), size = 2) +
        scale_y_log10() +
        scale_x_log10() +
        scott_theme_1() +
        annotation_logticks(sides = 'lb') +
        labs(x = 'p-value of top prediction', y = 'Number of conditions', title = 'Top gene set prediction\nper compound') +
        theme(axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
              axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
              plot.title = element_text(margin = margin(0, 0, 15, 0))
              )


    pdf(file = file.path(pre, sprintf('%s.pdf', filename_meat)), height = 6, width = 6)
    print(dvp_plot)
    dev.off()

    # Print the values
    write_vals = function(tab, sample_type, filename_meat, pre) {
        filename = file.path(pre, sprintf('%s_%s.txt', filename_meat, sample_type))
        write.table(tab, file = filename, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
    top_counts_dt[, write_vals(.SD, .BY[[1]], filename_meat, pre), by = sample_type]

}

# Define function to enforce monotonicity
enforce_monotonicity = function(vec) {

    # Since the first value(s) can be NA
    # (0 control / 0 treatment), start
    # after that.
    first_non_na_idx = which(!is.na(vec))[1]
    non_na_vec = vec[first_non_na_idx:length(vec)]
    for (i in 1:(length(non_na_vec) - 1)) {
        if (non_na_vec[i + 1] < non_na_vec[i]) {
            non_na_vec[i + 1] = non_na_vec[i]
        }
    }

    vec[first_non_na_idx:length(vec)] = non_na_vec
    vec
}

# Define function to adjust the FDR to the smallest observed
# FDR at any further point in the vector (inspired heavily
# by Benjamini-Hochberg FDR procedure)
smooth_fdr = function(vec) {

    # FDR values are provided in increasing order based on p-value,
    # so they are already "sorted" but not completely. To deal with
    # this, reverse the order and compute the cumulative minimum
    # value.
    rev_vec = rev(vec)
    rev_cummin = cummin(rev_vec)

    # There can be NAs at the end of this list that we want to
    # inherit the minimum value!
    rev_cummin[is.na(rev_cummin)] = min(rev_cummin, na.rm = TRUE)

    rev(rev_cummin)
}

get_pval_fdr_mapping <- function(top_dat) {

    # Here I get a table that says, "For every unique p-value, this is how many
    # predictions are made."
    get_top_pval_last_obs_dt <- function(samp_type, top_dat) {
        top_dat[sample_type == samp_type][, list(max_scaled_count = max(scaled_condition_count)), by = p_value]
    }

    join_pval_last_obs_dts_to_all_pvals <- function(top_dat_list, pvals) {
        for (i in 1:length(top_dat_list)) {
            setkey(top_dat_list[[i]], p_value)
            top_dat_list[[i]] <- top_dat_list[[i]][J(pvals), roll = TRUE]
            top_dat_list[[i]][is.na(max_scaled_count), max_scaled_count := 0]
        }

        top_dat_list
    }

    calc_fdr <- function(control, treatment, top_dat_list) {
        fdr <- top_dat_list[[control]][['max_scaled_count']] / top_dat_list[[treatment]][['max_scaled_count']]
        # If num discovered treatment conds is zero but it is above zero
        # for the control, set FDR to 1. This can be brought down later
		# as determined by the "smooth_fdr" function.
        fdr[(top_dat_list[[treatment]][['max_scaled_count']] == 0) & (top_dat_list[[control]][['max_scaled_count']] > 0)] <- 1
		# Chop FDR values above 1 to just 1.
        fdr[fdr > 1] <- 1

        print('making sure fdr is indeed sorted!')
        print('original:')
        print(fdr[1:50])
        print('fdr smoothed (now guaranteed monotonic)')
        fdr = smooth_fdr(fdr)
        print(fdr[1:50])
       
        # Commenting out since this should not be necessary.
        ## Since the FDR list is sorted here, I can also add an important step:
        ## not letting the FDR go down after it reaches 100%
        #first_100_fdr <- which(fdr == 1)[1]
        #fdr[first_100_fdr:length(fdr)] <- 1

        fdr
    }

    calc_fdr_all_controls <- function(top_dat_list, controls, treatment) {
        lapply(controls, calc_fdr, treatment = treatment, top_dat_list = top_dat_list)
    }

    make_fdr_only_dt <- function(fdrs_by_control, pvals) {

        get_single_fdr_dt <- function(fdr_type, fdrs_list, pvals) {
            data.table(p_value = pvals, fdr = fdrs_list[[fdr_type]], fdr_type = fdr_type)
        }

        fdr_dts <- lapply(names(fdrs_by_control), get_single_fdr_dt, fdrs_list = fdrs_by_control, pvals = pvals)

        rbindlist(fdr_dts)
    }

    
    ##################################
    #####   End subfunctions   #######
    #####   Begin "script"     #######
    ##################################
   
    setkeyv(top_dat, c('p_value'))
    
    sample_types <- unique(top_dat[['sample_type']])
    treatment <- 'treatment'
    controls <- sample_types[sample_types != treatment]

    # Here we get one table per sample type. Each table shows how many predictions are made
    # at each unique p-value.
    top_pval_last_obs_dt_list <- lapply(sample_types, get_top_pval_last_obs_dt, top_dat = top_dat)

    # Here all p-values from the tables immediately above are combined into a master pool
    # of p-values (and sorted, of course).
    all_pvals_unsorted <- Reduce(union, lapply(top_pval_last_obs_dt_list, `[[`, 'p_value'))
    all_pvals_sorted <- sort(all_pvals_unsorted)

    # This adds all p-values from the recently created p-value pool into each of the
    # "top_pval_last_obs" tables. This enables calculation of FDR for each control
    # in a consistent manner. The "roll" action of data.table allows the number of
    # counted conditions ("max_scaled_count") to be filled in appropriately when a
    # p-value is joined that is not currently in the table (it inherits from the row
    # above).
    top_pval_last_obs_dt_list <- join_pval_last_obs_dts_to_all_pvals(top_pval_last_obs_dt_list, all_pvals_sorted)
    names(top_pval_last_obs_dt_list) <- sample_types

    # The FDR estimate is actually calculated here. The control table and treatment table are
    # already lined up, so just divide the number of control conditions at a particular p-value
    # by the number of treatment conditions at the same p-value. Gold.
    print(top_pval_last_obs_dt_list)
    fdrs_by_control <- calc_fdr_all_controls(top_pval_last_obs_dt_list, controls = controls, treatment = treatment)
    names(fdrs_by_control) <- controls

    make_fdr_only_dt(fdrs_by_control, all_pvals_sorted)

}

add_fdr<- function(dat, pval_fdr_map, controls) {

    pval_fdr_map_wide <- data.table(dcast(pval_fdr_map, p_value ~ fdr_type, value.var = 'fdr'))
    new_control_type_names = sprintf('%s_FDR', controls)
    setnames(pval_fdr_map_wide, controls, new_control_type_names)
    
    setkey(dat, p_value)
    setkey(pval_fdr_map_wide, p_value)

    for (i in seq_along(new_control_type_names)) {
        control_type = new_control_type_names[i]
        dat = pval_fdr_map_wide[, c('p_value', control_type), with = FALSE][dat, roll = TRUE]
    }

    setkey(dat, p_value)

    dat
}

plot_fdr_vs_discoveries <- function(top_dt_with_fdr, controls, filename_meat, pre) {
   
    print(controls)

    split_off_last = function(x, delim) {
        split_x = strsplit(x, delim)[[1]]
        paste(split_x[-length(split_x)], collapse = delim)
    }

    stack_fdr_cols <- function(top_dt) {
    
        new_control_type_names = sprintf('%s_FDR', controls)
        dt_base <- top_dt[, !(names(top_dt) %in% new_control_type_names), with = FALSE]
        long_dt_list <- lapply(new_control_type_names, function(x) data.table(dt_base, fdr = top_dt[[x]], fdr_type = split_off_last(x, '_')))
        rbindlist(long_dt_list)
    }
    
    top_dt_fdr_long <- stack_fdr_cols(top_dt_with_fdr)
   
    # Do the plot thing
    fvd_plot <- ggplot(data = top_dt_fdr_long[sample_type == 'treatment']) +
        geom_step(aes(x = scaled_condition_count, y = fdr, color = fdr_type), size = 2, direction = 'vh') +
        scale_x_log10() +
        scott_theme_1() +
        labs(x = 'Number of conditions', y = 'False Discovery Rate', title = 'False discovery rate based on the\ntop prediction for each condition') +
        scale_colour_brewer(palette = 'Set1') +
        annotation_logticks(sides = 'b') +
        theme(axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
              axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
              plot.title = element_text(margin = margin(0, 0, 15, 0))
              )
    
    pdf(file.path(pre, sprintf('%s.pdf', filename_meat)), height = 6, width = 6)
    print(fvd_plot)
    dev.off()
    
    # Print the values
    write_vals = function(tab, fdr_type, filename_meat, pre) {
        print(fdr_type)
        filename = file.path(pre, sprintf('%s_%s.txt', filename_meat, fdr_type))
        write.table(tab, file = filename, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
    top_dt_fdr_long[sample_type == 'treatment', write_vals(.SD, .BY[[1]], filename_meat, pre), by = fdr_type]
}

write_pval_zscore_dt_old <- function(p_z_dt, filename_meat, pre) {

    p_z_dt <- p_z_dt[, .SD[order(p_value, -z_score)], by = condition]

    fname = file.path(pre, sprintf('%s.txt.gz', filename_meat))
    f = gzfile(fname, 'wb')
    write.table(p_z_dt, f, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
    close(f)
}

write_pval_zscore_dt <- function(p_z_dt, filename_meat, pre) {

    p_z_dt <- p_z_dt[order(p_value, -z_score)]

    fname = file.path(pre, sprintf('%s.txt.gz', filename_meat))
    f = gzfile(fname, 'wb')
    write.table(p_z_dt, f, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
    close(f)
}

write_top_pval_zscore_dt <- function(p_z_top_dt, filename_meat, pre) {

    p_z_top_dt <- p_z_top_dt[order(p_value, -z_score)]

    fname = file.path(pre, sprintf('%s.txt', filename_meat))
    f = file(fname, 'wt')
    write.table(p_z_top_dt, f, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
    close(f)
}

write_top_pval_zscore_dt_glob_by_sample_type <- function(p_z_top_dt, filename_meat, pre) {

    p_z_top_dt <- p_z_top_dt[, .SD[order(p_value, -z_score)], by = sample_type]

    con_name <- 'file'
    filename <- file.path(pre, sprintf('%s.txt', filename_meat))

    table_to_glob(p_z_top_dt, filename, by_vec = 'sample_type', connection_FUN = con_name)

}
