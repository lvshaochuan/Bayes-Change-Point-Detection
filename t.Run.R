  CP <- res$CP
  Q <- res$Q
  theta <- res$param
  Tau <- res$Tau
  Alpha <- res$Alpha
  Ui <- res$Ui
  changepoints <- Changepoints <- intensity <- jumprates <- CPnum <- Taus <- Alphas <- NULL
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

  for (j in 1:iter){
  if (res$Number != 1){  
  rates <- -diag(Q)
  rates <- rates[-length(rates)]
  rates <- c(rates, 1/(ending-starting))
  flag <- FALSE
  while (flag==FALSE){
  A <- Uniform(CP, rates)
  Uniformtimes <- A$U
  index <- which(diff(c(starting,Uniformtimes,ending)) <= diff(eventtimes[1:2]))
  Uniformtimes <- Uniformtimes[-index]
  JumpProb <- A$JumpProb[-index, ]
  if ((length(Uniformtimes) >=3) && (length(Uniformtimes) <= 300)) flag <- TRUE
  }
  }
  if (res$Number == 1){
  
  flag <- FALSE
  while (flag==FALSE){
  A <- Uniform(c(starting, ending), 1/(ending-starting))
  Uniformtimes <- A$U
  index <- which(diff(c(starting,Uniformtimes,ending)) <= diff(eventtimes[1:2]))
  Uniformtimes <- Uniformtimes[-index]
  JumpProb <- A$JumpProb[-index, ]
  if ((length(Uniformtimes) >=3) && (length(Uniformtimes) <= 300)) flag <- TRUE
  }
  }
  print(j)
 # FFBS.t <- function(Uniformtimes, JumpProb, eventtimes, obs, param, DF, Tau, Ui, Alpha)
  res <- FFBS.t(Uniformtimes, JumpProb, eventtimes, obs, param=theta, DF=1, Tau, Ui, Alpha)
  num <- res$Number
  cp <- res$cp
  CP <- res$CP
  Q <- res$Q 
  KK <- length(CP)
  if ((ending == CP[KK-1])&&(res$Number != 1)){
  CP <- CP[-KK]
  cp <- cp[-length(cp)]
  Q <- Q[,-(KK-1)]
  Q <- Q[-(KK-2),]
  num <- num - 1
  
  }  
  theta <- res$param
  Tau <- res$Tau
  Alpha <- res$Alpha
  Ui <- res$Ui
  if (j > burnin){
  changepoints <- c(changepoints, list(cp))
  Changepoints <- c(Changepoints, list(CP))
  jumprates <- c(jumprates, list(Q))
  CPnum <- c(CPnum, num)
  intensity <- c(intensity, list(theta))
  Taus <- c(Taus, Tau)
  Alphas <- c(Alphas, Alpha)
        }
  }      
    