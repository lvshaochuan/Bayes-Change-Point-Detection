  CP <- res$CP
  Q <- res$Q
  theta <- res$theta
  Var <- res$Var
  changepoints <- Changepoints <- intensity <- Variance <- jumprates <- CPnum <- NULL
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
  if ((max(diff(Uniformtimes) < 3.9)) && (length(Uniformtimes) <= 300)  && (length(Uniformtimes) >=3)) flag <- TRUE
#  if ((length(Uniformtimes) >=3) && (length(Uniformtimes) <= 300)) flag <- TRUE
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
  if ((max(diff(Uniformtimes) < 3.9)) && (length(Uniformtimes) <= 300)  && (length(Uniformtimes) >=3)) flag <- TRUE
#  if ((length(Uniformtimes) >=3) && (length(Uniformtimes) <= 300)) flag <- TRUE
  }
  }
  print(j)
  
  res <- FFBS.Objective(Uniformtimes, JumpProb, eventtimes, obs, theta, Var)
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
  theta <- res$theta
  Var <- res$Var
  if (j > burnin){
  changepoints <- c(changepoints, list(cp))
  Changepoints <- c(Changepoints, list(CP))
  jumprates <- c(jumprates, list(Q))
  CPnum <- c(CPnum, num)
  intensity <- c(intensity, list(theta))
  Variance <- c(Variance, Var)
        }
  }      
 
 
 ###############################################################################
 ###############################################################################
  position <- NULL
 for (i in 1:length(CPnum)){
 temp <- length(Changepoints[[i]])
 position <- c(position, Changepoints[[i]][-c(1,temp)])
 } 
 
  thetas <- NULL
 for (i in 1:length(CPnum)){
 temp <- length(intensity[[i]])
 thetas <- c(thetas, intensity[[i]][-c(1,temp)])
 }     