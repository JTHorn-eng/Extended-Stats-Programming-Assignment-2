# Task 1
# This is a final working version.
# This vector is the first 1000 elements of up to 5000 elements, where each value 1:n has an
# equal chance of being repeated 1:hmax times.
# Indices of the vector that have the same value share the same house.
n = 10000; hmax = 5
h1 <- rep(rep(1:n, sample(1:hmax,n, replace=TRUE)), length.out=n)

# Task 2
get.net = function(beta, h, n_c=15) {

    # Setup M
    M <- matrix(h,n,n)
    Mt <- t(M)
    M <- 1*(Mt == M)

    beta_mean <- mean(beta)

    # Compute the replacement matrix
    # Use outer product of beta vector
    sociability_matrix <- (n_c * outer(beta, beta)) / (beta_mean^2 * (n - 1))

    # Fill new matrix with runif values
    random_matrix <- matrix(0, n, n)
    random_matrix <- runif(n * n)

    # Compare the matrix of probabilities
    result_matrix <- matrix(0, n, n)
    result_matrix[sociability_matrix > random_matrix] <- 2

    # Filter upper triangle of M by the upper triangle of the results_matrix
    mask <- upper.tri(result_matrix) & result_matrix == 2
    M[mask] <- result_matrix[mask]

    # Mirror the upper tri to the lower
    M[lower.tri(M)] <- t(M)[lower.tri(M)]

    # Output a list of indices
    return(lapply(seq_len(nrow(M)), function(i) which(M[i, ] == 2)))

}


output <- get.net(runif(n), h1,15)