# Function to add a pre-filtered (to only negative controls) dummy
# dataset sample table to the original cg sample table for which
# predictions will be made.
add_dummy_controls_to_cg_sample_tab= function(cg_sample_tab, cg_control_col, cg_name_col, dummy_control_sample_tab, dummy_control_col, dummy_name_col) {

    dummy_control_sample_tab = copy(dummy_control_sample_tab)
    # Change names if necessary so the control columns line up
    if (!(is.null(cg_control_col) | is.null(dummy_control_col))) {
        setnames(dummy_control_sample_tab, dummy_control_col, cg_control_col)
    } else if (is.null(cg_control_col) & (!is.null(dummy_control_col))) {
        stop('A negative control column for the dummy dataset sample table was specified, but none specified for the original dataset. Please specify a negative control column for the original dataset sample table and try again.')
    } else if (is.null(dummy_control_col)) {
        return(cg_sample_tab)
    }
    if (!(is.null(cg_name_col) | is.null(dummy_name_col))) {
        setnames(dummy_control_sample_tab, dummy_name_col, cg_name_col)
    } else if (is.null(cg_name_col) & (!is.null(dummy_name_col))) {
        warning('A condition name column for the dummy dataset sample table was specified, but none was specified for the original dataset. It will not be used.')
    }

    rbindlist(list(cg_sample_tab, dummy_control_sample_tab), fill = TRUE, use.names = TRUE)

}

