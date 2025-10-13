

# ask simon, if you need to resample
###########################################################
# File:    Submission.R
# Authors: Colin McDonagh (s2213748), Matthew Poole (s2892092)
# , James Horn (s2304075)
# Date:    2025-10-13
###########################################################
#
### About
#
# This document contains an implementation of an SEIR model as 
# described in the Extended Statistical Programming module.
#
#
### Contributions:
#
# All members attended the majority of workshops and communicated remotely about progress made outside of
# class time.
#
# Colin McDonagh (%)  
#                      
#                     
#                      
#                      
#
# Matthew Poole  (%)  
#                      
#                      
#                      
#                       
#
# James Horn     (%)  
#                     
#                     
#

n = 10000; hmax = 5
h1 <- rep(rep(1:n, sample(1:hmax,n, replace=TRUE)), length.out=n)
beta <- runif(n)

get.net = function(beta, h, n_c=15) {
  
  # Setup M

  M <- outer(h1,h1,FUN="==")
  
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


nseir <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
  n <- length(h)
  ni <- round(n * pinf)

  k <- (alpha[3]*nc*beta)/(mean(beta)^2*(n-1))
  
  x <- rep(0, n)
  
  start_inf <- sample(1:n, ni, replace=FALSE)
  x[start_inf] <- 2
  
  S <- E <- I <- R <- t_vec <- rep(0,nt)
  S[1] <- sum(x==0)
  I[1] <- sum(x==2)

  pre_unif1 <- matrix(runif(n * nt), nrow = nt, ncol = n)
  pre_unif2 <- matrix(runif(n * nt), nrow = nt, ncol = n)
  
  for (t in 2:nt){

    any_inf_vec <- rep(0,n)
    
    ii_inf <- which(x==2)
    
    for (i in ii_inf){
      house_ii <- which(h==h[i])
      pnet_ii <- alink[[i]]
      
      personal_con <- rep(0,n)
      personal_con[house_ii] <- alpha[1]
      personal_con[pnet_ii] <- alpha[2]
      
      rnet_con <- beta[i] * k
      
      unif1 <- pre_unif1[, t]
      unif2 <- pre_unif2[, t]
      
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

plot_seir = function(seir_results, title = 'SEIR Model Results') {

    par(mfcol=c(1,1),mar=c(4,4,2,1))

    n_days <- length(seir_results$S)
    days <- 1:n_days

    max_y <- max(c(seir_results$S, seir_results$E, seir_results$I), na.rm = TRUE)

    plot(
        days
        , seir_results$S
        , type = "l"
        , col = "black"
        , ylim = c(0, max_y)
        , xlab = "Day"
        , ylab = "Number of People"
    )
    lines(days, seir_results$E, col = "blue")
    lines(days, seir_results$I, col = "red")
    lines(days, seir_results$R, col = "green")

    legend("topright", legend = c("S (Susceptible)", "E (Exposed)", "I (Infectious)", "R (Recovered)"),
            col = c("black", "blue", "red", "green "), lty = 1, cex = 0.8)

    mtext(title, outer = TRUE, line = -1.5, cex = 1.2)

}



# Run Simulation
alink <- get.net(runif(n), h1,15)
nseir(beta, h1, alink)
