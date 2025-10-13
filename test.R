# Task 1
# This is a final working version.
# This vector is the first 1000 elements of up to 5000 elements, where each value 1:n has an
# equal chance of being repeated 1:hmax times.
# Indices of the vector that have the same value share the same house.
n = 10000; hmax = 5
h1 <- rep(rep(1:n, sample(1:hmax,n, replace=TRUE)), length.out=n)

# Task 2
# I have started task 2. M outputs a matrix where aij = 1 if i and j share a household, 0 otherwise

print(h1)



n_c <- 15
beta <- runif(n)
beta_mean <- mean(beta)

# Apply sociability probabilities

# Using normal for looping for every matrix ij entry
# is very slow with n=10000, therefore:
#
# Only update 0 upper triangle values and use
# vectorised computation

# Identify zero entries only in the upper triangle (excluding diagonal)
# zero_upper <- (M == 0)

# Compute the replacement matrix
# Use outer product of beta vector
replacement <- (n_c * outer(beta, beta)) / (beta_mean^2 * (n - 1))
replacement
# Replace only those zeros in the upper triangle
# M[zero_upper] <- replacement[zero_upper]

# # Mirror upper triangle to lower
# M[lower.tri(M)] <- t(M)[lower.tri(M)]

# Return list of vector indexes   
# a <- lapply(seq_len(ncol(M)), function(i) M[,i])


