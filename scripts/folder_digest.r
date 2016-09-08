#!/usr/bin/env Rscript

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script is either run straight from the command line or is
# sourced into another R script. Its purpose is to take a folder,
# recursively generate a list md5 checksums for all of the
# contents of that folder, and then compute an md5 checksum on
# that list. This is for computing checksums on datasets present
# as multiple files in a folder, as comparing checksums of
# tarred, gzipped versions of those folders has proven difficult.
#
# For example, I have dataset in folder X and I want to check it
# against tarred & gzipped folder Y in the online repository.
# If I tar & gzip folder X to compare it to Y, there's no
# guarantee that the checksums will match, as I have discovered
# that R's internal tar function is the most reliable way to
# tar things from R, yet it is not convenient for creating new
# tarred folders to upload to the data repository. Therefore,
# the user (or me) should use this script to compute the md5
# checksum on the untarred, uncompressed contents of the data
# folder in question, and then deposit this value in the
# dataset configuration table. Comparisons will work
# seamlessly thereafter!

library(digest)

get_folder_md5 = function(folder) {
    
    files = list.files(folder, recursive = TRUE, include.dirs = FALSE, full.names = TRUE)
    md5_vec = vapply(files, function(x) digest(file = x, algo = 'md5'), character(1), USE.NAMES = FALSE)
    
    files_md5_formatted = paste(files, md5_vec, sep = '\t', collapse = '\n')
    message(files_md5_formatted)
    
    # Return md5 of the md5 vector!
    digest(md5_vec, algo = 'md5')
}

# If this is run from the command line instead of sourced from
# another script where the '__main__' option is set to FALSE,
# then this "main" function will run. Otherwise, this script
# is only good for the function it provides.
if(getOption('__main__', default = TRUE)) {

    folder = commandArgs(TRUE)
    #print(folder)
    final_md5 = get_folder_md5(folder)
    
    # Write to stdout!
    write(final_md5, file = '')

}
