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

get.net = function(beta, h, n_c=15) {
  #' Creates a list of vectors, with each vector at position i containing the indices of people that person i is connected with socially.
  #' This is done by looping through each person, calculating their network based on the probability they are connected to someone else, 
  #' and then recording these values at position i of the list, and then recording the value i at each index of the list of the person
  #' in the social network.
  #' @param beta A list of uniform random deviates
  #' @param h a vector describing who shares houses
  #' @param n_c a value dictating the average number of connections per person.
  #' 
  #' @return a list of vectors, with the vector at position i being the social network of person i.
  
  # Define the number of people to consider
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


nseir <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005) {
  #' This function implements and simulates the SEIR model described in the introduction to this document for a given set of parameter
  #' values. It calculates the amount of people that are in the susceptible, exposed, infected, and recovered classes each day based 
  #' on households, a personal social network, and a random interaction network (where households and the personal networks are disjoint). 
  #' @param beta A list of uniform random deviates
  #' @param h A vector indexing who shares houses
  #' @param alink A list of vectors, with the vector at position i being the social network of person i.
  #' @param alpha A vector of parameters for alpha_h, alpha_c, and alpha_r (as described above) respectively.
  #' @param delta The probability that an infected individual transfers from the infected class to the recovered class each day.
  #' @param gamma The probability that an exposed individual transfers from the exposed class to the infected class each day.
  #' @param n_c A value of the average number of personal connections per person.
  #' @param nt The number of days to simulate the model over.
  #' @param pinf The proportion of the population infected at the start of the simulation.
  #' 
  #' @return A list of vectors, where these vectors contain the number of people in the susceptible class, exposed class, infected class, 
  #' recovered class, and the day these correspond to respectively.
  
  # Sets the number of people in the model
  n <- length(h)
  
  # Finds the initial number of people who should be infected
  ni <- round(n * pinf)
  
  # Creates a vector to contain the class of each individual each day
  x <- rep(0, n)
  
  # Calculates a constant used inside a loop once for efficiency
  k <- (alpha[3]*nc*beta) / (mean(beta)^2*(n-1))
  
  # Randomly samples the initial infected people from our population
  start_inf <- sample(1:n, ni, replace=FALSE)
  x[start_inf] <- 2
  
  # Initialise the vectors we will return from this function
  S <- E <- I <- R <- t_vec <- rep(0,nt)
  
  # Find the initial number of susceptible and infected people
  S[1] <- sum(x==0)
  I[1] <- sum(x==2)
  
  # Iterate the following over each day in the simulation
  for (t in 2:nt){
    
    # Initialise the vector to store people infected each day
    any_inf_vec <- rep(0,n)
    
    # Find the number of people infected on day t
    ii_inf <- which(x==2)
    
    # Iterate over each infected person to find out which susceptible people are infected by at least one infected person
    for (i in ii_inf){
      
      # Find the household and personal social network links of the current infected person
      house_ii <- which(h==h[i])
      pnet_ii <- alink[[i]]
      
      # Create a new vector where the jth entry will be the probability that person i infects person j by household or by the personal network
      personal_con <- rep(0,n)
      
      # If person j is in the household network of person i, set the jth value in personal_con to alpha_h
      personal_con[house_ii] <- alpha[1]
      # If person j is in the household network of person i, set the jth value in personal_con to alpha_c
      personal_con[pnet_ii] <- alpha[2]
      
      # Find the probabilities of person i infecting each other person in the random network.
      rnet_con <- k * beta[i]
      
      # Find which people are conncected to person i by household or personal network
      con_ii <- which(personal_con != 0)
      
      # Make a vector of uniform random deviates where only indices that correspond to contacts are sampled
      unif1 <- rep(1,n)
      unif1[con_ii] <- runif(length(con_ii))
      
      # Separate uniform random deviate of length n
      unif2 <- runif(n)
      
      # Simulate the probability that each person was infected by person i by household or personal network
      inf_vec1 <- personal_con >= unif1
      # Simulate the probability that each person was infected by person i by random network
      inf_vec2 <- rnet_con >= unif2
      
      # Stores whether each person was infected by any infected person in the model on day t
      any_inf_vec <- inf_vec1 | inf_vec2 | any_inf_vec
    }
    
    # Create uniform random deviates
    u <- runif(n)
    
    # Infected people transition to recovered with class probability delta            
    x[x==2&u<delta] <- 3 
    # Exposed people transition to infected class with probability gamma
    x[x==1&u<gamma] <- 2  
    
    # (Only) Susceptible people transition to exposed class as calculated above
    x[x==0&(1*any_inf_vec)==1] <- 1
    
    # Sum the total people in each class and increment the day
    S[t] <- sum(x==0); E[t] <- sum(x==1)
    I[t] <- sum(x==2); R[t] <- sum(x==3)
    t_vec[t] <- t
  }  
  
  # Return the list of susceptible, exposed, infected, and recovered people on each day and the vector of day indices
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
  
  # Create subplot titles
  plot_titles <- c(
    'Standard Model'
    ,'Random Mixing Only'
    ,'Constant Sociability Parameters'
    ,'Random Mixing Only and Constant Sociability Parameters'
  )
  
  for (i in seq_along(seir_results_list)) {
    
    # Obtain results timeframe in days
    seir_results <- seir_results_list[[i]]
    n_days <- length(seir_results$S)
    days <- 1:n_days
    
    # Generate line graphs for each class member of the population per day
    plot(
      days
      , seir_results$S
      , type="l"
      , col = "black"
      , ylim = c(0, n)
      , xlab = "Day"
      , ylab = "Number of People"
      , main = plot_titles[i]
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
  mtext(title, outer = TRUE, line = 1, cex = 1.2)
}


# Run models
# Specify the number of people in the simulation and the maximum number of people in a given household
system.time({n = 10000; hmax = 5

# Simulate beta values as standard uniform random deviates
beta <- runif(n)

# Simulate a vector of household contacts
h1 <- rep(rep(1:n, sample(1:hmax,n, replace=TRUE)), length.out=n)

# Simulate personal network information
alink <- get.net(runif(n), h1, 15)

# Create a list of models (using nseir) to compare
results_list <- list(
  
  # Standard model with personal, household, and random mixing
  nseir(beta, h1, alink)
  
  # Run model with random mixing only
  ,nseir(beta, h1, alink, alpha=c(0, 0, 0.04))
  
  # Run model with constant sociability parameters for each person
  ,nseir(beta=rep(mean(beta),n), h1, alink)
  
  # Run model with constant sociability parameters for each person and only consider random mixing
  ,nseir(beta=rep(mean(beta),n), h1, alink, alpha=c(0, 0, 0.04))
  
)
plot_seir(results_list)
})
# Standard Model
# -------------------------------------------------------------------------

# This represents a typical SEIR model with household and personal 
# networks. Day 0-10 a small number of infections initially appear. Between
# days 30-40 the number of infected people reaches its peak as many susceptible
# people are already infected/immune. As a result, between day 50-70 
# new infections decrease rapidly and finally the epidemic stablises around day 
# 80-100 as almost all people have recovered.


# Random mixing only
# -------------------------------------------------------------------------

# We see a similar trend compared with the standard model
# in terms of changes in SEIR numbers, however the number of exposed 
# people peaks sooner and larger at around day 20.

# Infected individuals mostly contact people they already know
# within household/personal networks, meaning once their household 
# or close network becomes infected, they run out of new susceptibles
# before the infection reaches the rest of the population.

# However without personal and household networks, the population 
# mixes homogenously with less specific population interaction
# ,causing the E and I curves to rise and peak earlier. It appears
# that the rate of recovery is also larger, however this is due to larger
# numbers of infected people per day, the rate of recovery remains constant
# in both this model and the standard.
# 
