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