# Functions are defined first.
# See below functions to see the model
set.seed(0)
get.net = function(beta, h, nc=15) {
  
  # Setup M
  M <- matrix(h,n,n)
  Mt <- t(M)
  M <- 1*(Mt == M)
  
  beta_mean <- mean(beta)
  
  # Compute the replacement matrix
  # Use outer product of beta vector
  sociability_matrix <- (nc * outer(beta, beta)) / (beta_mean^2 * (n - 1))
  
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

nseir1 <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
  n <- length(h)
  ni <- round(n * pinf)
  
  x <- rep(0, n)
  
  start_inf <- sample(1:n, ni, replace=FALSE)
  x[start_inf] <- 2
  
  S <- E <- I <- R <- t_vec <- rep(0,nt)
  S[1] <- sum(x==0)
  I[1] <- sum(x==2)
  
  # constant vector to calculate random infected probabilities quicker
  # Don't know if this actually makes it quicker
  k <- (0.01*nc*beta)/(mean(beta)*(n-1))
  
  for (t in 2:nt){
    u <- runif(n) ## probabilities
    
    # indices of people living with infected person
    hi <- which(h %in% h[which(x==2)])
    
    # indices of people living in infected by household and get exposed
    hie <- hi[x[hi]==0 & u[hi]<alpha[1]] # replace 0.1 with alpha[1]
    
    # all indices of x in network of infecteds
    # what about people with multiple infected people in their network?
    ni <- unlist(alink[which(x==2)])
    
    # indices of people with infected person in their network and get exposed
    nie <- ni[x[ni]==0 & u[ni]<alpha[2]]
    
    # probability of susceptible person being exposed by a random infected person
    rprobs <- outer(beta[which(x==2)],k[which(x==0)],FUN="*")
    
    # realising these probabilities and finding indices of people exposed
    # what if we get integer(0) because noone gets infected randomly? Think it's fine as it gives number(0)
    rie <- which(colSums(rprobs > runif(length(rprobs)))>0)

x[x==2&u<delta] <- 3 ## I -> R with prob delta 
x[x==1&u<gamma] <- 2 ## E -> I with prob gamma 
x[c(hie,nie,rie)] <- 1 ## S -> E 
S[t] <- sum(x==0); E[t] <- sum(x==1)
I[t] <- sum(x==2); R[t] <- sum(x==3)
t_vec[t] <- t
  }  
  return(list(S=S,E=E,I=I,R=R,t=t_vec))
}
#####################################################################
# Start by generating houses
# This vector is the first 1000 elements of up to 5000 elements, where each value 1:n has an
# equal chance of being repeated 1:hmax times.
n = 10000; hmax = 5
house_info <- rep(rep(1:n, sample(1:hmax,n, replace=TRUE)), length.out=n)

# generate sociabilities that should remain constant throughout I think
sociabilities = runif(n)

# create the network using get.net()
# This is still taking a long time
network_list <- get.net(beta=sociabilities,h=house_info)

#run the model using nseir1()
# This is now very speedy :)
output <- nseir1(beta=sociabilities,h=house_info,alink=network_list)
output

