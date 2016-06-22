set_prefix <- function(pre) {
    dir.create(pre, recursive = TRUE)
    pre
}

file_timestamp <- function() {

    now <- as.POSIXlt(Sys.time())
    
    with(now, sprintf('%04d-%02d-%02d_%02d.%02d', 1900 + year, mon + 1, mday, hour, min))
}

