# Task 1
# This is a final working version.
# This vector is the first 1000 elements of up to 5000 elements, where each value 1:n has an
# equal chance of being repeated 1:hmax times.
# Indices of the vector that have the same value share the same house.
n = 10000; hmax = 5
h1 <- rep(rep(1:n, sample(1:hmax,n, replace=TRUE)), length.out=n)

# Task 2
# I have started task 2. M outputs a matrix where aij = 1 if i and j share a household, 0 otherwise

get.net = function(beta, h, n_c=15) {

    M <- matrix(h,n,n)
    Mt <- t(M)
    M <- 1*(Mt == M)

    beta_mean <- mean(beta)

    # Apply sociability probabilities

    # Using normal for looping for every matrix ij entry
    # is very slow with n=10000, therefore:
    #
    # Only update 0 upper triangle values and use
    # vectorised computation

    # Identify zero entries only in the upper triangle (excluding diagonal)
    zero_upper <- (M == 0) & upper.tri(M)

    # Compute the replacement matrix
    # Use outer product of beta vector
    replacement <- (n_c * outer(beta, beta)) / (beta_mean^2 * (n - 1))

    # Replace only those zeros in the upper triangle
    M[zero_upper] <- replacement[zero_upper]

    # # Mirror upper triangle to lower
    M[lower.tri(M)] <- t(M)[lower.tri(M)]

    # Return list of vector indexes    
    return(lapply(seq_len(ncol(M)), function(i) M[,i]))

}

# Takes up to 10 seconds
task2_ls <- get.net(runif(n), h1,15)


nseir = function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005) {


    print('[INFO]  Running SEIR model for t days')
    
    S <- E <- I <- R <- rep(0,n)

    # For each day ...
    for (t in 1:nt) {









    }


}

nseir(
    beta = beta
    ,h = h1
    ,alink = get.net(beta, h1)
)