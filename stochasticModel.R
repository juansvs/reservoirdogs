library(tidyverse)

####### functions #####
# Rates for 6 possible events: transmission to sp.1, transmission to sp.2
# recovery of sp.1, recovery of sp.2, shedding, mortality
nxtvnt <- function(X) {
  # give names to the vector entries
  S1 <- X[1]; S2 <- X[2]
  I1 <- X[3]; I2 <- X[4]
  E <- X[5]
  # calculate the mean rate at which each event is happening
  trans1 <- b1*S1*E
  trans2 <- b2*S2*E
  recov1 <- r1*I1
  recov2 <- r2*I2
  shed <- l1*I1+l2*I2
  mrtl <- E*g
  # combine the rates into a single vector
  rates <- c(trans1, trans2, recov1, recov2, shed, mrtl)
  # calculate the times to the next events;
  # draw 6 random numbers from a uniform distribution, determine the log
  # and multiply by the inverse of the corresponding rate
  deltas <- -1/rates*log(runif(6))
  return (list(event=which.min(deltas), time=min(deltas)))
}


# Create object to store results
res <- data.frame(beta2 = betas[1], 
                  iter=1, 
                  t = 0, 
                  S1 = n0[1], 
                  S2 = n0[2], 
                  I1 = n0[3], 
                  I2 = n0[4], 
                  E = n0[5])

fadeoutT <- data.frame(beta2 = 0, t = 0)

##### changing values of beta ####
microbenchmark::microbenchmark(
  for (b2 in betas) {
    for (i in 1:nsims) {
      if (i>1) {
        res <- rbind(res, c(i, 0, n0))
      }
      
      # Set initial objects that will be updated as the simulation progresses
      N <- n0
      t <- 0
      
      counter <- 0
      day <- 1
      while (t<730) {
        event <- nxtvnt(N)
        
        updtVec <- switch (event$event,
                           c(-1, 0, 1, 0, 0),
                           c(0, -1, 0, 1, 0),
                           c(1, 0,-1, 0, 0),
                           c(0, 1, 0, -1, 0),
                           c(0, 0, 0, 0, 1),
                           c(0, 0, 0, 0, -1)
        )
        
        N <- N+updtVec
        if (N[5]==0) {fadeoutT <- rbind(fadeoutT, c(b2, t)); break}
        t <- t+event$time
        # if (event$event <=4) {
        #   res <- rbind(res, c(i, t, N))
        # }
        counter <- counter + 1
        if (ceiling(t)>day) {
          day <- ceiling(t)
          res <- rbind(res, c(b2, i, t, N))
          
        }
        # day <- ceiling(t)
        # if (counter %% 50 ==0) {
        #   cat("Iteration:",i, "Step:", counter, "time =", t, "Pop sizes =",N, "\n")
        # }
      }
    }
    
  }
  )

##### different initial populations sizes ####
resN0 <- data.frame(n2_0 = N2s[1], iter=1, t = 0, S1 = n0[1], S2 = N2s[1], I1 = n0[3], I2 = n0[4], E = n0[5])
fadeoutN0 <- data.frame(n0 = N2s[1], t = 0)
b2 <- betas[2]
# parameter loop. Change some parameter and do several iterations with that parameter value

for (no in N2s) {
  n0 <- c(S1 = 29, S2 = no, I1 = 1, I2 = 0, E = 0)
  # replicates loop, do the same thing for a number of stochastic simulations
  for (i in 1:nsims) {
    if (i>1) {
      resN0 <- rbind(resN0, c(i, 0, n0))
    }
    
    # Set initial objects that will be updated as the simulation progresses
    N <- n0
    t <- 0
    
    counter <- 0
    day <- 1
    # time loop, while time is below two years, keep calculating the next event and updating the population vector accordingly
    while (t<730) {
      event <- nxtvnt(N)
      
      updtVec <- switch (event$event,
                         c(-1, 0, 1, 0, 0),
                         c(0, -1, 0, 1, 0),
                         c(1, 0,-1, 0, 0),
                         c(0, 1, 0, -1, 0),
                         c(0, 0, 0, 0, 1),
                         c(0, 0, 0, 0, -1)
      )
      
      N <- N+updtVec
      if (N[5]==0) {fadeoutN0 <- rbind(fadeoutN0, c(no, t)); break}
      t <- t+event$time
      # if (event$event <=4) {
      #   res <- rbind(res, c(i, t, N))
      # }
      counter <- counter + 1
      if (ceiling(t)>day) {
        day <- ceiling(t)
        resN0 <- rbind(resN0, c(no, i, t, N))
        
      }
      # day <- ceiling(t)
      # if (counter %% 50 ==0) {
      #   cat("Iteration:",i, "Step:", counter, "time =", t, "Pop sizes =",N, "\n")
      # }
    }
  }
  
}

######### analysis #######
res %>% group_by(iter) %>% filter(t>700) %>% transmute(prev1 = I1/(S1+I1), prev2 = I2/(S2+I2)) %>% ungroup() %>% 
  summarize(p1 = mean(prev1), sd1 = sd(prev1), p2 = mean(prev2), sd2 = sd(prev2))
# quantify the proportion of disease fading out.
res %>% filter(t == max(t)) %>% group_by(iter) %>% select(t, I1, I2) 

# calculate the proportion of fadeouts, and the mean time to extinction of the parasite
fadeoutT %>% group_by(beta2) %>% summarise(meant=mean(t), sdt=sd(t), N = length(t), prop = length(t)/nsims)

