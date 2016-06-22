# This script generates a function to compute cosine similarity between
# each column in one matrix with each column in another matrix


cos_sim <- function(x, y, na.rm = FALSE) {

    # Calculates cosine similarity between the columns of two matrices

    # If x or y are vectors, coerce to one-column matrix
    x <- as.matrix(x)
    y <- as.matrix(y)

    # Stop if x or y has one row - cannot compute similarities!
    stopifnot(dim(x)[1] > 1)
    stopifnot(dim(y)[1] > 1)
    
    # Set missing values in each to 0 (if na.rm == TRUE, for speed)
    if (na.rm) {
        x[is.na(x)] <- 0
        y[is.na(y)] <- 0
    }

    # Normalize x and y matrices (by column, of course)
    x_norms <- apply(x, 2, function(z) norm(as.matrix(z), 'f'))
    y_norms <- apply(y, 2, function(z) norm(as.matrix(z), 'f'))

    normed_x <- sweep(x, 2, x_norms, '/')
    normed_y <- sweep(y, 2, y_norms, '/')

    # Calculate cosine similarity and return the similarity matrix!
    t(normed_x) %*% normed_y

}
