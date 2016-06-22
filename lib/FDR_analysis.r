library(data.table)
library(Cairo)
library(ggplot2)
library(reshape2)

source('/project/csbio/Scott/Software/myRLibrary/plotting_tools/scott_themes.r')
source('/project/csbio/Scott/Software/myRLibrary/Munging/table_to_globular.r')


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
    
    dvp_plot <- ggplot(data = top_counts_dt) +
        geom_line(aes(x = p_value, y = scaled_condition_count, color = sample_type), size = 2) +
        scale_y_log10() +
        scale_x_log10() +
        scott_theme_1() +
        annotation_logticks(sides = 'lb') +
        labs(x = 'p-value of top prediction', y = 'Number of conditions', title = 'Top gene set prediction p values')


    CairoPDF(file = file.path(pre, sprintf('%s.pdf', filename_meat)), height = 6, width = 6)
    print(dvp_plot)
    dev.off()

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
        fdr[fdr > 1] <- 1

        print('making sure fdr is indeed sorted!')
        print(fdr[1:50])

        # Since the FDR list is sorted here, I can also add an important step:
        # not letting the FDR go down after it reaches 100%
        first_100_fdr <- which(fdr == 1)[1]
        fdr[first_100_fdr:length(fdr)] <- 1

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
    
    stack_fdr_cols <- function(top_dt) {
    
        new_control_type_names = sprintf('%s_FDR', controls)
        dt_base <- top_dt[, !(names(top_dt) %in% new_control_type_names), with = FALSE]
        long_dt_list <- lapply(new_control_type_names, function(x) data.table(dt_base, fdr = top_dt[[x]], fdr_type = strsplit(x, '_')[[1]][1]))
        rbindlist(long_dt_list)
    }
    
    top_dt_fdr_wide <- stack_fdr_cols(top_dt_with_fdr)
    
    fvd_plot <- ggplot(data = top_dt_fdr_wide[sample_type == 'treatment']) +
        geom_step(aes(x = scaled_condition_count, y = fdr, color = fdr_type), size = 2, direction = 'vh') +
        scale_x_log10() +
        scott_theme_1() +
        labs(x = 'Number of discoveries', y = 'False Discovery Rate', title = sprintf('FDR as assessed by\n%s conditions', paste(controls, collapse = ' and '))) +
        scale_colour_brewer(palette = 'Set1') +
        annotation_logticks(sides = 'b')
    
    CairoPDF(file.path(pre, sprintf('%s.pdf', filename_meat)), height = 6, width = 6)
    print(fvd_plot)
    dev.off()
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
