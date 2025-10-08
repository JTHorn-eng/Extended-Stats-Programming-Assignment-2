# Task 1
# This is a final working version.
# This vector is the first 1000 elements of up to 5000 elements, where each value 1:n has an
# equal chance of being repeated 1:hmax times.
# Indices of the vector that have the same value share the same house.
n = 10000; hmax = 5
h1 <- rep(rep(1:n, sample(1:hmax,n, replace=TRUE)), length.out=n)

# Task 2
# I have started task 2. M outputs a matrix where aij = 1 if i and j share a household, 0 otherwise
M <- matrix(h1,n,n)
Mt <- t(M)
M <- 1*(Mt == M)



