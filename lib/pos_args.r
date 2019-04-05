#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

library(optparse)

format_positional_arguments = function(pos_args_list) {

    arg_lengths = nchar(names(pos_args_list))
    longest_arg_len = max(arg_lengths)
    hanging_indent = longest_arg_len + 4
    spacers = vapply(arg_lengths, function(x) paste(character(hanging_indent - x + 1), collapse = ' '), character(1))
    combined_string = paste(names(pos_args_list), '%s', pos_args_list, '\n\n', sep = '', collapse = '\n')
    wrapped_string = strwrap(combined_string, width = 50 + hanging_indent, exdent = hanging_indent, prefix = '  ')
    wrapped_string = paste(wrapped_string, collapse = '\n')
    wrapped_with_spacers = do.call(sprintf, args = c(list(fmt = wrapped_string), as.list(spacers)))
    return(wrapped_with_spacers)
}

get_usage_and_positional = function(pos_args_list) {
    
    formatted_positional_arguments = format_positional_arguments(positional_arguments_list)

    usage_string = paste('%prog [options] ', paste(names(positional_arguments_list), collapse = ' '), sep = '')
    positional_arguments_help = paste('positional arguments:', formatted_positional_arguments, sep = '\n')

    final_usage_positional = sprintf('%s\n\n%s', usage_string, positional_arguments_help)
    return(final_usage_positional)

}
