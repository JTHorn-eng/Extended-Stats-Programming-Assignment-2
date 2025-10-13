nseir <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
  n <- length(h)
  ni <- round(n * pinf)
  
  x <- rep(0, n)
  
  start_inf <- sample(1:n, ni, replace=FALSE)
  x[start_inf] <- 2
  
  S <- E <- I <- R <- t_vec <- rep(0,nt)
  S[1] <- sum(x==0)
  I[1] <- sum(x==2)
  
  for (t in 2:nt){

    any_inf_vec <- rep(0,n)
    
    ii_inf <- which(x==2)
    
    for (i in ii_inf){
      house_ii <- which(h==h[i])
      pnet_ii <- alink[[i]]
      
      personal_con <- rep(0,n)
      personal_con[house_ii] <- alpha[1]
      personal_con[pnet_ii] <- alpha[2]
      
      rnet_con <- (alpha[3]*nc*beta*beta[i])/(mean(beta)^2*(n-1))
      
      unif1 <- runif(n)
      unif2 <- runif(n)
      
      inf_vec1 <- personal_con >= unif1
      inf_vec2 <- rnet_con >= unif2
      
      inf_vec <- inf_vec1 | inf_vec2
      
      any_inf_vec <- inf_vec | any_inf_vec
    }
    
    u <- runif(n) ## uniform random deviates
    x[x==2&u<delta] <- 3 ## I -> R with prob delta 
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma 
    x[x==0&(1*any_inf_vec)==1] <- 1 ## S -> E with prob beta*I[i-1] 
    S[t] <- sum(x==0); E[t] <- sum(x==1)
    I[t] <- sum(x==2); R[t] <- sum(x==3)
    t_vec[t] <- t
  }  
  return(list(S=S,E=E,I=I,R=R,t=t_vec))
}



# Task 1
# This is a final working version.
# This vector is the first 1000 elements of up to 5000 elements, where each value 1:n has an
# equal chance of being repeated 1:hmax times.
# Indices of the vector that have the same value share the same house.
n = 10000; hmax = 5
h1 <- rep(rep(1:n, sample(1:hmax,n, replace=TRUE)), length.out=n)

# Task 2
# I have started task 2. M outputs a matrix where aij = 1 if i and j share a household, 0 otherwise

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


alink <- get.net(runif(n), h1,15)


beta <- runif(n)
nseir(beta, h1, alink)