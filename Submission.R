###########################################################
# File:    Submission.R
# Authors: Colin McDonagh (s2213748), Matthew Poole (s2892092), James Horn (s2304075)
# Date:    2025-10-15
# Github URL: https://github.com/JTHorn-eng/Extended-Stats-Programming-Assignment-2/blob/main/Submission.R
###########################################################
#
### About
#
# This document is the code to run an SEIR model and prints the plots of 4 
# different implementations of the model. The program primarily uses two 
# functions to do this: get.net() and seir(). This SEIR model can account for 
# different transmission rates between people sharing houses of different sizes,
# social networks of different sizes, and random mixing. The housing sharing 
# assumes a uniform distribution of house sizes and the social networks allow 
# for variable social network sizes, depending on their 'sociability' parameter.
# The code concludes with a commentary on the 4 different models. 
# 
#
#
### Contributions:
#
# All members attended the majority of workshops and communicated remotely about 
# progress made outside of class time.
#
# Colin McDonagh (33 %) Wrote the majority of the nseir function and led the 
#                       commentary of the different models  
#                      
# Matthew Poole  (34 %) Wrote the code to distribute people to houses, optimised 
#                       the get.net function and nseir function
#                      
# James Horn     (33 %) Wrote the original get.net function and wrote the plot 
#                       function 
#                     


get.net = function(beta, h, nc=15) {
  #' Creates a list of vectors, with each vector at position i containing the 
  #' indices of people that person i is connected with socially. This is done by 
  #' looping through each person, calculating their network based on the 
  #' probability they are connected to someone else, and then recording these 
  #' values at position i of the list, and then recording the value i at each 
  #' index of the list of the person in the social network.
  #' @param beta A list of n uniform random deviates
  #' @param h    A vector of length n describing who shares houses
  #' @param nc  A value dictating the average number of connections per person.
  #' 
  #' @return a list of vectors, with the vector at position i being the social 
  #' network of person i.
  
  # Define the number of people to consider
  n <- length(h)
  
  # Create a list of vectors for who shares a house with who, as people sharing 
  # houses cannot be in the same social network.
  house_network <- lapply(h, function(x) which(x == h))
  
  # Generate an empty list to populate with connections
  social_network <- vector("list",n)
  
  # Calculates a constant used inside a loop once for efficiency to calculate 
  # the probability that two people are connected socially.
  k <- nc*beta/(mean(beta)^2*(n-1))
  
  # loop through each person in the network
  for (i in 1:(n-1)) {
    # Create probabilities that person i is connected with each other person. 
    # Only consider person i and onwards as these values will be symmetric.
    network_probs <- c(rep(0,i), beta[i]*k[(i+1):n]) 
    
    # Create probabilities to compare the network probs with.
    # The first i values are 1 as the first (i-1) have already been assigned 
    # networks, and person i cannot be in a network with themself. 
    probs <- c(rep(1,i), runif(n-i))
    
    # Assign the indices of people sharing a network to a vector, and exclude 
    # anyone who lives with them
    social_indices <- which(network_probs > probs)[-house_network[[i]]]
    
    # append the indices of people in person i's network to the list
    social_network[[i]] <- c(social_network[[i]], social_indices)
    
    # append person i to the network of each person in person i's network
    social_network[social_indices] <- lapply(social_network[social_indices], 
                                             function(x) c(x, i))
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
  
  # Initialise the vectors we will return from this function
  S <- E <- I <- R <- t_vec <- rep(0,nt)
  
  # Creates a vector to contain the class of each individual each day
  x <- rep(0, n)
  
  # Randomly samples the initial infected people from our population
  x[sample(1:n, round(n * pinf))] <- 2
  
  # Find the initial number of susceptible and infected people
  S[1] <- sum(x==0); I[1] <- sum(x==2)
  
  # Calculates a constant used inside a loop once for efficiency
  k <- (alpha[3]*nc*beta) / (mean(beta)^2*(n-1))
  
  # is it necessary to do this with k or is it worth for performance?
  
  # Iterate the following over each day in the simulation
  for (t in 2:nt){
    
    # Initialise a vector to store the people exposed by each infected person
    exp_vector <- vector("list", length=length(which(x==2)))
    
    # Iterate over each infected person to find out which susceptible people are
    # exposed by at least one infected person
    for (i in which(x==2)){
      
      # Reset the indices of susceptible people infected after each iteration
      shce <- snce <- srce <- numeric(0)
      
      # Only consider the household model if the user specifies this
      if(alpha[1]>0){
        # Find indices of susceptible people in a household with infected person i
        shc <- which((h==h[i])&x==0)
        # Samples whether people in shc are exposed (by household) by person i
        shce <- shc[runif(length(shc))<=alpha[1]]
      }
      
      
      # Only consider the contact network if the user specifies this
      if(alpha[2]>0){
        # Find indices of susceptible people in a contact network with infected 
        # person i
        nc <- alink[[i]]
        snc <- nc[which(x[nc]==0)]
        
        # Samples whether people in snce are exposed (by contact network) by person i
        snce <- snc[runif(length(snc))<=alpha[2]]
      }
      
      
      # Only consider the random network if the user specifies this
      if(alpha[3]>0){
        # Finds indices of susceptible people
        src <- which(x==0)
        # Samples whether people in src are exposed (by random network) by person i
        srce <- src[runif(length(src)) <= beta[i]*k[src]]
      } 
      
      # Append unique people exposed by person i to a running list
      exp_vector[[i]] <- unique(c(shce, snce, srce))
    }
    
    # Find the unique susceptible people exposed on day t by any infected person
    any_exp_vector <- unique(unlist(exp_vector))
    
    # Create uniform random deviates
    u <- runif(n)
    
    # Infected people transition to recovered with class probability delta            
    x[x==2&u<delta] <- 3 
    # Exposed people transition to infected class with probability gamma
    x[x==1&u<gamma] <- 2  
    
    # Susceptible people transition to exposed class as calculated above
    x[any_exp_vector] <- 1
    
    # Sum the total people in each class and increment the day
    S[t] <- sum(x==0); E[t] <- sum(x==1)
    I[t] <- sum(x==2); R[t] <- sum(x==3)
    t_vec[t] <- t
  }  
  
  # Return the list of susceptible, exposed, infected, and recovered people on each day and the vector of day indices
  return(list(S=S,E=E,I=I,R=R,t=t_vec))
}

plot_seir = function(seir_results_list, n, title = 'SEIR Model Results') {

    #' Plots a list of SEIR model results in a 2x2 grid layout, showing the
    #' evolution of SEIR (Susceptible, Exposed, Infectious, Recovered) 
    #' compartments for one or more simulation results. Each element of
    #' `seir_results_list` should be a list containing vectors `S`, `E`,
    #' `I`, and `R` representing daily counts of individuals in each state.
    #' 
    #' @param seir_results_list A list of all model results
    #' @param n the population size for each plot's y-axis
    #' @param title optional title for overall window

    # Generate the grid layout with margin outside layout for legend
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

    # Create subplot titles
    plot_titles <- c(
      'Standard Model'
      ,'Random Mixing Only'
      ,'Constant Sociability Parameters'
      ,'Random Mixing Only and Constant Sociability Parameters'
    )

    for (i in seq_along(seir_results_list)) {

        par(xpd=FALSE)
        
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
            , main = plot_titles[i]
        )

        grid()
      
        lines(days, seir_results$E, col = "blue")
        lines(days, seir_results$I, col = "red")
        lines(days, seir_results$R, col = "green")

        legend("topright",
              inset = c(0, 0),
              legend = c("S (Susceptible)", "E (Exposed)", "I (Infectious)", 
                         "R (Recovered)"),
              col = c("black", "blue", "red", "green"),
              lty = 1,
              xpd = NA,
              bty = "n",
              cex = 0.8
          )
    }

    # End layout plotting
    par(xpd = NA)

    # Add a title for the overall window
    mtext(title, outer = TRUE, line = 1, cex = 1.2)
}

system.time({
# Run models
# Specify the number of people in the simulation and the maximum number of 
# people in a given household
n = 10000; hmax = 5

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
    
    # Run model with constant sociability parameters for each person and only 
    # consider random mixing
    ,nseir(beta=rep(mean(beta),n), h1, alink, alpha=c(0, 0, 0.04))

)
plot_seir(results_list, n)

})

################################################################################
#'
###' Model Comparison
#'
##' Household, Contact Network, and Random Network Models vs Random Network Only Models:
# ------------------------------------------------------------------------------
#' The models that contain the household and contact network components account
#' for more variability in the infection rates. Thus, we observe a more diffuse 
#' distribution of infectious individuals than the models that only consider the 
#' random network effect, which overestimate the peak severity of the pandemic. 
#' We also observe that the peak of infection occurs later for models that contain
#' household and contact network components than it does for the random network
#' models. 
#' 
#' The additional variability described above comes from allowing the models that 
#' contain household and personal network information to consider behavioral and 
#' random aspects of individuals more intimately:
#' 
#' Variability is introduced in the personal network since the more social an 
#' individual is, the more people they are likely to have in their contact 
#' network, and thus the more people they are likely to infect/be infected by 
#' (where as the opposite is true for unsociable people).
#' 
#' The household structure adds variability that is independent of sociability
#' parameter values (unlike the random network only model for variable beta 
#' values) as household sizes are assumed to be uniformly distributed from 1 to 5.
#' 
#' We also note that the full model is more robust when neglecting the random 
#' variability of sociability between individuals, as the results are very similar
#' (only around 100 more people are infected in the constant sociability model).
#' This could be the case as the household and contact network structures still 
#' encode some underlying random sociability structure into the model). 
#' The random network only model is far more sensitive to this change, and 
#' results in the pandemic infecting far more people (and being severe for 
#' longer) than the random network model with variable sociability parameters.
# ------------------------------------------------------------------------------
