library(tidyverse)
library(doParallel)
####### functions #####
# Rates for 4 possible events: transmission, recovery, appearance of feces in the environment, and disappearance of feces

nxtvnt <- function(X) {
  # give names to the vector entries
  S <- X[1] 
  I <- X[2]
  E <- X[3]
  # calculate the mean rate at which each event is happening
  trans <- beta*S*E
  recov <- rho*I
  shed <- lambda0
  mrtl <- E*mu
  # combine the rates into a single vector
  rates <- c(trans, recov, shed, mrtl)
  # calculate the times to the next events;
  # draw 4 random numbers from a uniform distribution, determine the log
  # and multiply by the inverse of the corresponding rate
  deltas <- -1/rates*log(runif(4))
  return (list(event=which.min(deltas), time=min(deltas)))
}

updtVec <- function(event) {
  switch (event$event,
          c(-1, 1, 0), # transmission
          c(1, -1, 0), # recovery, adult parasite mortality
          c(0, 0, 1), # scat deposition
          c(0, 0, -1)# scat disappearance/loss of infectiousness
          )
}
  

#### Load parameters ####
source('params2.R')
beta_orig <- beta
mu_orig <- mu
mb_scenarios <- expand.grid(mu=mu_orig*deltas, beta=beta_orig*deltas)[-1,]

#### setup parallel computation ####
cl <- makeCluster(6)
registerDoParallel(cl)

#### Simulation output ####
# Open output file to store results. We can use a summarized output that only retains the first time to invasion
outfile <- paste(Sys.Date(),'_','simOutput.csv',sep='')
# outInvT <- 'stochastic_output_time.csv'
outfileconn <- file(outfile, open = 'w'); cat('beta','mu','iter', 'prev','tinv', file = outfileconn, sep=',', fill = T)
# outInvTconn <- file(outInvT, open = 'w'); cat('beta','iter','t', file = outInvTconn, sep = ',', fill = T)

#### run simulation ####
# bencmarking
# ptime <- system.time({
  # to run in parallel use this line:
  # res <- foreach (d = seq_len(length.out = dim(mb_scenarios)[1]), .combine = rbind) %dopar% {
for (d in seq_len(length.out = dim(mb_scenarios)[1])) {
  beta <- mb_scenarios[d,'beta']
  mu <- mb_scenarios[d,'mu']
  # foreach (i = 1:nsims, .combine = cbind) %dopar% {
  for (i in 1:nsims) {
    # Set initial objects that will be updated as the simulation progresses
    N <- n0
    t <- 0
    # cat(beta,i, t, N, file = outfileconn, sep=',', fill = T) # write initial state to file
    counter <- 0
    day <- 1
    invaded <- F 
    tinv <- integer()
    while (t<sims.duration) {
      event <- nxtvnt(N)
      N <- N+updtVec(event)
      t <- t+event$time
      if (event$event==1 && !invaded) {
        tinv <- t
        invaded <- T
      }

      # counter <- counter + 1
      # if (ceiling(t)>day) {
      #   day <- ceiling(t)
      #   cat(beta,i, round(t,2),N,file = outfileconn, sep=',', fill = T)
      # }
    }
    prev <- N[2]/(N[1]+N[2])
    cat(beta, mu, i, round(prev, 2), round(tinv, 2), file=outfileconn, sep=',', fill = T)
    # c(beta, mu, i, round(prev, 2), round(tinv, 2))#for parallel version
  }
};close(outfileconn)


# })[3]
ptime
stopCluster(cl)

#### analyze output ####
# stochastic prevalence mean as a function of beta*lambda0 and mu*rho
read_csv(outfile, col_types = 'ddddd') %>% 
  group_by(beta, mu) %>% summarise(meanprev=mean(prev),meant=mean(tinv, na.rm = T)) %>% 
  ggplot(aes(beta*lambda0/(beta_orig*lambda0), mu*rho/(mu_orig*rho)))+
  geom_contour_filled(aes(z=meanprev))+
  theme_minimal(base_size = 14)+
  labs(fill='Mean Prevalence', x=expression(paste(beta,lambda[0])),y=expression(paste(mu,rho)))

# time to invasion as a function of beta*lambda0 and mu*rho 
g3 <- read_csv(outfile, col_types = 'ddddd') %>% 
  group_by(beta, mu) %>% summarise(meanprev=mean(prev),meant=mean(tinv, na.rm = T)) %>% 
  ggplot(aes(beta*lambda0/(beta_orig*lambda0), mu*rho/(mu_orig*rho)))+
  geom_contour_filled(aes(z=meant))+
  theme_minimal(base_size = 12)+
  labs(fill='Mean time', x=expression(paste(beta,lambda[0])),y=expression(paste(mu,rho)))

# 95% CI of invasion times
g4 <- read_csv(outfile, col_types = 'ddddd') %>% 
  group_by(beta, mu) %>% summarise(mint=quantile(tinv, 0.025), maxt=quantile(tinv,0.975),tdif=maxt-mint) %>% 
  ggplot(aes(beta*lambda0/(beta_orig*lambda0), mu*rho/(mu_orig*rho)))+
  geom_contour_filled(aes(z=tdif))+
  theme_minimal(base_size = 12)+
  labs(fill='Time range', x=expression(paste(beta,lambda[0])),y=expression(paste(mu,rho)))

# the only case where there is no invasion is if beta =0, otherwise, 100% of cases result in invasion.
g2 <- read_csv(outfile, col_types = 'ddddd') %>% 
  group_by(beta, mu) %>% summarise(meanprev=mean(prev), rangeprev=quantile(prev,0.975)-quantile(prev, 0.0225)) %>% 
  ggplot(aes(beta*lambda0/(beta_orig*lambda0), mu*rho/(mu_orig*rho)))+
  geom_contour_filled(aes(z=rangeprev), binwidth = 0.1)+
  theme_minimal(base_size = 12)+
  # theme(aspect.ratio = 1)+
  labs(fill='Prev. range', x=expression(paste(beta,lambda[0])),y=expression(paste(mu,rho)))
g2
# range of prevalence
read_csv(outfile, col_types = 'ddddd') %>% 
  mutate(bdif=beta-beta_orig,mdif=mu-mu_orig) %>% 
  filter(bdif==(abs(bdif)), mdif==min(abs(mdif))) %>% summarise(range(prev))

# proportion of invasions
read_csv(outfile, col_types = 'ddddd') %>% 
  group_by(beta, mu) %>% mutate(inv=tinv>0) %>% 
  mutate(bdif=beta-beta_orig,mdif=mu-mu_orig) %>% 
  filter(bdif==min(abs(bdif)), mdif==min(abs(mdif))) %>% summarise(range(prev))
