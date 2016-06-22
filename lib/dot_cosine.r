# Here I define Raamesh's 'dotcosine' method for comparing chemical genetic to
# genetic interaction profiles.  It normalizes the columns of the second matrix
# but not those of the first.


dot_cos_sim <- function(x, y, na.rm = FALSE) {

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

    # Normalize only y matrix! (by column, of course)
    y_norms <- apply(y, 2, function(z) norm(as.matrix(z), 'f'))

    normed_y <- sweep(y, 2, y_norms, '/')

    # Calculate dot-cosine similarity and return the similarity matrix!
    t(x) %*% normed_y

}
