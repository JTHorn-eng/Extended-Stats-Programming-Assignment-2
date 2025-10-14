system.time({
# Task 1
# This is a final working version.
# This vector is the first 1000 elements of up to 5000 elements, where each value 1:n has an
# equal chance of being repeated 1:hmax times.
# Indices of the vector that have the same value share the same house.
n = 10000; hmax = 5
h1 <- rep(rep(1:n, sample(1:hmax,n, replace=TRUE)), length.out=n)


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
    network_probs <- c(rep(0,i-1), (beta[i]*nc*beta[i:n])/(mean(beta)^2*(n-1))) 
    
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
  
  social_network
  
}



task2_ls <- get.net(runif(n), h1,15)
})