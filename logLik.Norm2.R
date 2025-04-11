logLik.Norm2 <- function(eventtimes, obs, sim_mu, sim_sigma, starts, ends){
N <- 0
z <- 0
i <- 1
logLiks <- 0
###########  
while(eventtimes[i] <= starts) i <- i+1
########### <
while (eventtimes[i] < ends){
logLiks <- logLiks + dnorm(obs[i], mean=sim_mu, sd=sim_sigma, log=T)
  i <- i+1
  if (i > length(eventtimes)) break
  }
z <- logLiks
return(z)
}