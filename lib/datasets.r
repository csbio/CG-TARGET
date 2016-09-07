#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script contains functions that retrieve the filenames
# of the requested genetic interaction dataset and, if
# necessary, download said datasets from the specified
# location.
library(R.utils)
library(digest)

# Source in the function that generates an md5 checksum for
# a directory structure.
TARGET_PATH = Sys.getenv('TARGET_PATH')
options('__main__' = FALSE)
source(file.path(TARGET_PATH, 'scripts/folder_digest.r'))
options('__main__' = NULL)


url_username_password = function(url, username, password) {
    url_no_http = sub('^http://', '', url)
    
    sprintf('http://%s:%s@%s', username, password, url_no_http)
}

dir_file_preview = function(folder) {
    # Function to get the first 80 characters of the
    # first ten lines of each file in the folder for
    # previewing!
    files = list.files(folder, recursive = TRUE, include.dirs = FALSE, full.names = TRUE)

    f_preview_vec = character(length(files))
    for (i in 1:length(files)) {

        f = files[i]
        f_ten_lines = scan(f, what = character(), n = 10, sep = '\n')
        f_preview_vec[i] = paste(substr(f_ten_lines, 1, 80), collapse = '\n')

    }

    final_preview = paste(files, f_preview_vec, sep = '\n\n', collapse = '\n---------------------\n\n')

    # Return a string that contains one file preview
    # per file in the directory, pre-formatted with
    # newlines and everything!
    return(final_preview)

}



get_gi_info = function(gi_name, software_dir) {
   
    gi_data_dir = file.path(software_dir, 'data', 'GI')
    dataset_dir = file.path(gi_data_dir, gi_name)
    
    dataset_tab_filename = file.path(gi_data_dir, 'dataset_config.txt')
    dataset_tab = fread(dataset_tab_filename, header = TRUE, sep = '\t')
    setkey(dataset_tab, dataset_name)
    
    table_filename = file.path(dataset_dir, dataset_tab[gi_name][['matrix']])
    array_filename = file.path(dataset_dir, dataset_tab[gi_name][['row_info_table']])
    query_filename = file.path(dataset_dir, dataset_tab[gi_name][['col_info_table']])

    download_dataset = TRUE
    expected_md5 = dataset_tab[gi_name][['md5_checksum']]
    if (dir.exists(dataset_dir)) {
        downloaded_md5 = get_folder_md5(dataset_dir)

        print(sprintf('expected md5: %s', expected_md5))
        print(sprintf('downloaded md5: %s', downloaded_md5))

        if (downloaded_md5 == expected_md5) {
            # If the expected md5 matches with the one for the dataset
            # already downloaded, then there's nothing to do!
            download_dataset = FALSE
        }
    }

    if (download_dataset) {
        if (dir.exists(dataset_dir)) {
            # If the directory is there but checksums don't match up, then
            # remove the directory and start over!
            unlink(dataset_dir, recursive = TRUE)
        }
        
        download_gi_dataset(gi_name, gi_data_dir, dataset_tab, dataset_tab_filename)

        # Now check again to see if md5 matches! If not, abort!!!
        downloaded_md5 = get_folder_md5(dataset_dir)
        if (expected_md5 != downloaded_md5) {
            # Get the first 80 characters of the first ten lines of each file for previewing!
            preview = dir_file_preview(dataset_dir)
            #message(preview, '')

            unlink(dataset_dir)
            message('\n\n')
            download_location = dataset_tab[gi_name][['download_location']]
            stop(sprintf('\n\nUnable to download expected genetic interaction dataset from location:\n%s.\nThe md5 checksums of the expected and downloaded datasets do not match.\nPlease check that the username and password are correct and try again.\n\nPreviews of the offending file(s):\n\n%s\n\n', download_location, preview))
        }
    }

    #if (!all(file.exists(c(table_filename, array_filename, query_filename)))) {
    #    # If the data files are not present, then get them!
    #    # Gotta remove any existing directory with that name so the
    #    # new tarfile can be unzipped (R's internal method does not
    #    # appear to overwrite very well).
    #    if(dir.exists(dataset_dir)) {
    #        unlink(dataset_dir, recursive = TRUE)
    #        #file.remove(normalizePath(list.files(dataset_dir, full.names = TRUE)))
    #        #file.remove(dataset_dir)
    #    }
    #    download_gi_dataset(gi_name, gi_data_dir, dataset_tab, dataset_tab_filename)
    #    return(list(gi_tab = table_filename, gi_array_tab = array_filename, gi_query_tab = query_filename))
    #}

    # Get information for important dataset columns
    array_sys_name_col = dataset_tab[gi_name][['row_info_tab_sys_name_col']]
    query_genename_col = dataset_tab[gi_name][['col_info_tab_genename_col']]
    if(query_genename_col == '') {
        query_genename_col = NULL
    }
    
    return(list(gi_tab = table_filename, gi_array_tab = array_filename, gi_query_tab = query_filename,
                array_sys_name_col = array_sys_name_col, query_genename_col = query_genename_col))
}

download_gi_dataset = function(gi_name, gi_data_dir, dataset_tab, dataset_tab_filename) {

    download_location = dataset_tab[gi_name][['download_location']]
    if(is.null(download_location)) {
        stop(sprintf('Dataset "%s" not in directory "%s" and no download location was given in "%s"', gi_name, dataset_dir, dataset_tab_filename))
    }
    dest_file = file.path(gi_data_dir, 'tmp.tar.gz')
    
    # Ask for username/password
    message('About to download dataset at:\n', download_location)
    response = ''
    username = NULL
    password = NULL
    while(!(response %in% c('y', 'n'))) {
        message('Does the file require a username and password? (y/n)')
        response_con = file('stdin')
        response = readLines(response_con, 1)
        close(response_con)
        #response = scan(what = character())
    }
    if(response == 'y') {
        message('username:')
        username_con = file('stdin')
        username = readLines(username_con, 1)
        close(username_con)
        message('password:')
        password_con = file('stdin')
        password = readLines(username_con, 1)
        close(password_con)
        download_loc_userpwd = url_username_password(download_location, username, password)
    } else {
        download_loc_userpwd = download_location
    }

    download.file(download_loc_userpwd, destfile = dest_file)

    # downloadFile doesn't work for some reason now that I have
    # proper password authentication
    #downloadFile(url = download_location, filename = dest_file, username = username, password = password)

    ####### DO NOW
    # Provide some way to handle if the downloaded file cannot be untarred
    # (remove the temporary file, etc)
    untar_res = try(untar(dest_file, exdir = gi_data_dir, tar = 'internal'))
    if (class(untar_res) == 'try-error') {
        file_preview = scan(dest_file, n = 10, what = character(), sep = '\n')
        unlink(dest_file)
        message('\n\n')
        stop(sprintf('\n\nUnable to untar/decompress expected genetic interaction dataset from location:\n%s.\nPlease check that the username and password are correct and try again.\n\nPreview of offending genetic interaction tarred folder:\n\n%s\n\n', download_location, paste(file_preview, collapse = '\n')))
    }
    #gunzip(dest_file)
    unlink(dest_file)

}


get_gene_set_info = function(gi_name, gene_set_name, software_dir) {
   
    gene_set_dir = file.path(software_dir, 'data', 'gene_sets')
    gene_set_gi_dir = file.path(gene_set_dir, gi_name, gene_set_name)
    dir.create(gene_set_gi_dir, recursive = TRUE)
    
    gene_set_tab_filename = file.path(gene_set_dir, 'gene_set_config.txt')
    gene_set_tab = fread(gene_set_tab_filename, header = TRUE, sep = '\t')
    setkeyv(gene_set_tab, c('dataset_name', 'gene_set_name'))

    gene_set_gi_filename = file.path(gene_set_gi_dir, gene_set_tab[J(gi_name, gene_set_name)][['gene_set_filename']])
    gene_set_interpretable_column = gene_set_tab[J(gi_name, gene_set_name)][['interpretable_column']]
    if(is.na(gene_set_interpretable_column)) {
        gene_set_interpretable_column = NULL
    }
    
    # To determine if the gene-set file in question needs to 
    # be downloaded, first check to see if the file exists,
    # and then check to see if the checksum matches the
    # expected one.
    download_dataset = TRUE
    expected_md5 = gene_set_tab[J(gi_name, gene_set_name)][['md5_checksum']]
    if (file.exists(gene_set_gi_filename)) {
        downloaded_md5 = digest(file = gene_set_gi_filename, algo = 'md5')
        #print(expected_md5)
        #print(downloaded_md5)
        if (expected_md5 == downloaded_md5) {
            # Nothing to do if everything matches up :)
            download_dataset = FALSE
        }
    }
    if (download_dataset) {
        # Regardless of if the file doesn't exist or if it
        # does but doesn't match the checksum, download
        # it!!! This function will overwrite the file if
        # it exists
        download_location = gene_set_tab[J(gi_name, gene_set_name)][['download_location']]
        download_gene_set_table(gi_name, gene_set_name, gene_set_gi_dir, gene_set_gi_filename, download_location)
        
        # Now check again to see if md5 matches! If not, abort!!!
        downloaded_md5 = digest(file = gene_set_gi_filename, algo = 'md5')
        if (expected_md5 != downloaded_md5) {
            file_preview = scan(gene_set_gi_filename, n = 10, what = character(), sep = '\n')
            unlink(gene_set_gi_filename)
            message('\n\n')
            stop(sprintf('\n\nUnable to download expected gene-set table from location:\n%s.\nPlease check that the username and password are correct and try again.\n\nPreview of offending gene-set file:\n\n%s\n\n', download_location, paste(file_preview, collapse = '\n')))
        }
    }
    
    message(gene_set_gi_filename)
    message(gene_set_interpretable_column)
    
    return(list(filename = gene_set_gi_filename, interpretable_column = gene_set_interpretable_column))
}

download_gene_set_table = function(gi_name, gene_set_name, gene_set_gi_dir, gene_set_gi_filename, download_location) {

    #print(gene_set_tab)
    #print(gene_set_tab[J(gi_name, gene_set_name)])

    #download_location = gene_set_tab[J(gi_name, gene_set_name)][['download_location']]
    if(is.null(download_location)) {
        message('\n\n')
        stop(sprintf('File "%s" for gene set "%s" does not exist, and no download location was given in "%s"', gene_set_gi_filename, gene_set_name, gene_set_tab_filename))
    }
    
    # Ask for username/password
    message('About to download gene_set_file at:\n', download_location)
    response = ''
    username = NULL
    password = NULL
    while(!(response %in% c('y', 'n'))) {
        message('Does the file require a username and password? (y/n)')
        response_con = file('stdin')
        response = readLines(response_con, 1)
        close(response_con)
    }
    if(response == 'y') {
        message('username:')
        username_con = file('stdin')
        username = readLines(username_con, 1)
        close(username_con)
        message('password:')
        password_con = file('stdin')
        password = readLines(username_con, 1)
        close(password_con)
        download_loc_userpwd = url_username_password(download_location, username, password)
    } else {
        download_loc_userpwd = download_location
    }
    
    download.file(download_loc_userpwd, destfile = gene_set_gi_filename)
    
    #downloadFile(url = download_location, filename = gene_set_gi_filename, username = username, password = password, skip = FALSE)
    
    
    #untar(dest_file, exdir = gi_data_dir, tar = 'internal')
    #gunzip(dest_file)
    #file.remove(dest_file)

}
