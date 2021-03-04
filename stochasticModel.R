library(tidyverse)
library(deSolve)

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

ODEfunc <- function(t, states, params) {
  with(as.list(c(states, params)), {
    # rates of change
    dS1 <- -b1*S1*E+r1*I1
    dS2 <- -b2*S2*E+r2*I2
    dI1 <- b1*S1*E-r1*I1
    dI2 <- b2*S2*E-r2*I2
    dE <- l1*I1+l2*I2 - g*E
    
    list(c(dS1, dS2, dI1, dI2, dE))
  }
  )
}

# Parameter values
betas <- c(1e-7, 1e-6, 0.5e-5, 1e-5, 0.5e-4, 1e-4)
N2s <- c(50, 30, 20, 10)
b1 <- 0.1e-4; # b2 <- 0.1e-5
r1 <- 1/365; r2 <- 1/365
l1 <- 1; l2 <- 1
g <- 1/15

# Initial conditions. Vector in order: S1, S2, I1, I2, E
n0 <- c(S1 = 29, S2 = 20, I1 = 1, I2 = 0, E = 0)

# set number of iterations for each parameter combination
nsims <- 100


#### Deterministic model #######
for (b in 1:length(betas)) {
  b2 <- betas[b]
  # Specify parameters and initial conditions
  ODEparams <- c(b1, b2, r1, r2, l1, l2, g)
  ODEstates <- n0
  
  # specify times 
  times <- 0:730
  
  out <- ode(y = ODEstates, times = times, func = ODEfunc, parms = ODEparams)
  out <- as.data.frame(out)
  out$b2 <- b2
  if (b == 1) {
    outs <- out
  } else {    outs <- rbind(outs, out)}

}
###### Stochastic model ########

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
system.time(
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


######## plots ########
ggplot(res)+geom_line(aes(t, I1, group = iter), color = 'darkred', alpha = 0.4) +
  geom_line(aes(t, I2, group = iter), color='steelblue', alpha = 0.4)+
  geom_line(aes(time, I1), data = as_tibble(out), color = 'darkred',size = 1.5)+
  geom_line(aes(time, I2), data = as_tibble(out), color = 'steelblue', size = 1.5)+
  theme_classic(base_size = 16) +
  labs(x = 'Time (days)', y = 'Abundance')

# plot time to parasite extinction, proportion of extinctions, as a function of Beta2
fadeoutT %>% filter(beta2 != 0) %>% ggplot(aes(beta2/r2, t, group=beta2)) +geom_boxplot()+
  theme_minimal(base_size = 16)+
  labs(x = expression(paste(beta,2,'/',rho)), y = 'Time to extinction')
fadeoutT %>% filter(beta2 != 0) %>% group_by(beta2) %>% summarise(prop = length(t)/nsims) %>% 
  ggplot(aes(beta2/r2, prop)) + geom_point() + theme_classic(base_size = 16) +
  ylim(0,NA)+
  labs(x = expression(paste(beta,2,'/',rho)), y = 'Probability of extinction')

# plot time to parasite extinction, proportion of extinctions, as a function of N2_0
fadeoutN0 %>% ggplot(aes(n0, t, group=n0)) +geom_boxplot()+
  theme_minimal(base_size = 16)+
  labs(x = 'N2_0', y = 'Time to extinction')
fadeoutN0 %>% group_by(n0) %>% summarise(prop = length(t)/nsims) %>% 
  ggplot(aes(n0/30, prop)) + geom_point(size=2.5) + theme_classic(base_size = 16) +
  ylim(0,NA)+
  labs(x = 'N2_0/N1_0', y = 'Probability of extinction')


# Deterministic model prevalence in sp2 v time
outs %>% ggplot(aes(time, I2/(S2+I2)))+geom_point(aes(color = as.factor(b2/r2))) + theme_minimal(base_size = 16)+
  labs(x='Time(days)', y = "Prevalence", color=expression(paste(beta,2,'/',rho)))

# plot example trajectories with 2 different values of beta
res %>% filter(beta2 %in% betas[c(4,6)]) %>% 
  ggplot(aes(t, I2/(I2+S2), group=iter,color=as.factor(beta2)))+geom_path() +
  theme_classic(base_size = 16) +
  labs(x='Time (days)', y = 'Prevalence', color = expression(paste(beta,2)))#+
  geom_line(aes(time, I2/(S2+I2), color = as.factor(b2)), data=outs[outs$b2 %in% betas[c(4,6)],])

# select 10 iterations, example plot with 2 different values of beta2  
betaSamps <- betas[c(3,6)]
stochSamp1 <- res %>% filter(iter<=10, beta2 == betaSamps[1])
stochSamp2 <- res %>% filter(iter<=10, beta2 == betaSamps[2])
detSamp1 <- outs %>% filter(b2==betaSamps[1])
detSamp2 <- outs %>% filter(b2==betaSamps[2])
ggplot()+
  geom_line(aes(t, E, group = iter), data = stochSamp1, color = 'darkmagenta', alpha = 0.4)+
  geom_line(aes(t, E, group = iter), data = stochSamp2, color = 'darkgreen', alpha = 0.4)+
  geom_line(aes(time, E), data = detSamp1, color = 'darkmagenta', size = 1.5)+
  geom_line(aes(time, E), data = detSamp2, color = 'darkgreen', size = 1.5)+
  theme_minimal(base_size = 16)+
  labs(x = 'Time (days)', y = 'Abundance')

# select 10 iterations, example plot with 2 different values of N2_0  
n0samps <- N2s[c(1,4)]
stochSamp1 <- res %>% filter(iter<=10, n == betaSamps[1])
stochSamp2 <- res %>% filter(iter<=10, beta2 == betaSamps[2])
detSamp1 <- outs %>% filter(b2==betaSamps[1])
detSamp2 <- outs %>% filter(b2==betaSamps[2])
ggplot()+
  geom_line(aes(t, E, group = iter), data = stochSamp1, color = 'darkmagenta', alpha = 0.4)+
  geom_line(aes(t, E, group = iter), data = stochSamp2, color = 'darkgreen', alpha = 0.4)+
  geom_line(aes(time, E), data = detSamp1, color = 'darkmagenta', size = 1.5)+
  geom_line(aes(time, E), data = detSamp2, color = 'darkgreen', size = 1.5)+
  theme_minimal(base_size = 16)+
  labs(x = 'Time (days)', y = 'Abundance')

out %>% as_tibble() %>% ggplot(size = 2)+geom_line(aes(time, I1/(I1+S1))) +
  geom_line(aes(time, I2/(S2+I2)))+
  theme_classic(base_size = 16) +
  labs(x='Time (days)', y = 'Prevalence')

ggplot(res)+geom_line(aes(t, I1/(S1+I1), group = iter), color = 'darkred', alpha = 0.4) +
  geom_line(aes(t, I2/(S2+I2), group = iter), color='steelblue', alpha = 0.4)+
  geom_line(aes(time, I1/(S1+I1)), data = as_tibble(out), color = 'darkred',size = 1.5)+
  geom_line(aes(time, I2/(S2+I2)), data = as_tibble(out), color = 'steelblue', size = 1.5)+
  theme_classic(base_size = 16) +
  labs(x = 'Time (days)', y = 'Prevalence')


######### analysis #######
res %>% group_by(iter) %>% filter(t>700) %>% transmute(prev1 = I1/(S1+I1), prev2 = I2/(S2+I2)) %>% ungroup() %>% 
  summarize(p1 = mean(prev1), sd1 = sd(prev1), p2 = mean(prev2), sd2 = sd(prev2))
# quantify the proportion of disease fading out.
res %>% filter(t == max(t)) %>% group_by(iter) %>% select(t, I1, I2) 

# calculate the proportion of fadeouts, and the mean time to extinction of the parasite
fadeoutT %>% group_by(beta2) %>% summarise(meant=mean(t), sdt=sd(t), N = length(t), prop = length(t)/nsims)

