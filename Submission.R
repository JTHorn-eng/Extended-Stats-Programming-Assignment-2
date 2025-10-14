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
  #' Creates a list of vectors, with each vector at position i containing the indices of people that person i is connected with socially.
  #' This is done by looping through each person, calculating their network based on the probability they are connected to someone else, 
  #' and then recording these values at position i of the list, and then recording the value i at each index of the list of the person
  #' in the social network
  #' @param beta A list of uniform random deviates
  #' @param h a vector describing who shares houses
  #' @param n_c a value dictating the average number of connections per person.
  #' 
  #' @return a list of vectors, with the vector at position i being the social network of person i.
  
  
  n <- length(h)
  
  # Create a list of vectors for who shares a house with who, as people sharing houses cannot be in the same social network.
  house_network <- lapply(h, function(x) which(x == h))
  
  # Generate an empty list to populate with connections
  social_network <- vector("list",n)
  
  # loop through each person in the network
  for (i in 1:n) {
    # Create probabilities that person i is connected with each other person. 
    # Only consider person i and onwards as these values will be symmetric.
    network_probs <- c(rep(0,i-1), (beta[i]*n_c*beta[i:n])/(mean(beta)^2*(n-1))) 
    
    # Create probabilities to compare the network probs with.
    # The first i values are 1 as the first (i-1) have already been assigned networks, and person i cannot be in a network with themself. 
    probs <- c(rep(1,i), runif(n-i))
    
    # Assign the indices of people sharing a network to a vector, and exclude anyone who lives with them
    social_indices <- which(network_probs > probs)[-house_network[[i]]]
    
    # append the indices of people in person i's network to the list
    social_network[[i]] <- c(social_network[[i]], social_indices)
    
    # append person i to the network of each person in person i's network
    social_network[social_indices] <- lapply(social_network[social_indices], function(x) c(x, i))
  }
  
  return(social_network)
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

plot_seir = function(seir_results_list, title = 'SEIR Model Results') {

    #' Plots a list of SEIR model results in a 2x2 grid layout, showing the
    #' evolution of SEIR (Susceptible, Exposed, Infectious, Recovered) 
    #' compartments for one or more simulation results. Each element of
    #' `seir_results_list` should be a list containing vectors `S`, `E`,
    #' `I`, and `R` representing daily counts of individuals in each state.
    #' 
    #' @param seir_results_list A list of all model results
    #' @param title optional title for overall window
    #'     

    # Generate the grid layout with margin outside layout for legend
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 4, 6))

    for (i in seq_along(seir_results_list)) {
        
        # Obtain results timeframe in days
        seir_results <- seir_results_list[[i]]
        n_days <- length(seir_results$S)
        days <- 1:n_days

        # Generate line graphs for each class member of the population per day
        plot(
            days
            , seir_results$S
            , type = "l"
            , col = "black"
            , ylim = c(0, n)
            , xlab = "Day"
            , ylab = "Number of People"
        )
        lines(days, seir_results$E, col = "blue")
        lines(days, seir_results$I, col = "red")
        lines(days, seir_results$R, col = "green")
    }

    # End layout plotting
    par(xpd = NA)

    # Create a clear legend for all plots
    legend("topright",
        inset = c(-0.55, 0),
        legend = c("S (Susceptible)", "E (Exposed)", "I (Infectious)", "R (Recovered)"),
        col = c("black", "blue", "red", "green"),
        lty = 1,
        xpd = NA,
        bty = "n",
        cex = 0.8)

    # Add a title for the overall window
    mtext(title, outer = TRUE, line = -1.5, cex = 1.2)
}

# Run Simulation
alink <- get.net(runif(n), h1, 15)
results_list <- list(
    nseir(beta, h1, alink)

)
plot_seir(results_list)
