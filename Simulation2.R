 data1 <- rnorm(2048, mean=0, sd=10)
 data1[205:266] <- data1[205:266] + 14.64
 data1[267:307] <- data1[267:307] -  3.66
 data1[308:471] <- data1[308:471] + 7.32
 data1[472:511] <- data1[472:511] - 7.32
 data1[512:819] <- data1[512:819] + 10.98
 data1[820:901] <- data1[820:901] - 4.39
 data1[902:1331] <- data1[902:1331] + 3.29
 data1[1332:1556] <- data1[1332:1556] + 19.03
 data1[1557:1597] <- data1[1557:1597] + 7.68
 data1[1598:1658] <- data1[1598:1658] + 15.37
################################################################################

################################################################################
 savpar <- par(mfrow=c(3,1))
 plot(data1, cex=0.5, main="(a) Simulated Data", xlab="", ylab="")
 segments(1,0,204,0, col='red', lwd=2)
 segments(205,14.64,266,14.64, col='red', lwd=2)
 segments(267,-3.66,307,-3.66, col='red', lwd=2)
 segments(308,7.32,471,7.32, col='red', lwd=2)
 segments(472,- 7.32,511,- 7.32, col='red', lwd=2)
 segments(512,10.98,819,10.98, col='red', lwd=2)
 segments(820,- 4.39,901,- 4.39, col='red', lwd=2)
 segments(902,3.29,1331,3.29, col='red', lwd=2)
 segments(1332,19.03,1556,19.03,col='red', lwd=2)
 segments(1557,7.68,1597,7.68,col='red', lwd=2)
 segments(1598,15.37,1658,15.37, col='red', lwd=2)
 segments(1659,0,2048,0, col='red', lwd=2)
hist(position/eventtimes[1], 800, xlim=c(0, 2048), xlab="Changepoint Location", main="(b) Changepoint Locations")
box()
hist(CPnum[1001:length(CPnum)]-1, 80, xlab="Number of Changepoints", main="(c) Number of Changepoints", col="grey")
#hist(parameter, 500, xlab=expression(mu), main="(c) Histogram")
box()
par(savpar)
 
 
 source("Uniform.Adaptive.R")
 source("MarginallogLik.Norm.R")
 source("FFBS.Norm.R")

  starting <- 0
  ending <- 20.48
  eventtimes <- seq(0, ending, by=0.01) 
  eventtimes <- eventtimes[-1]
  obs <- (data1-5.9)/11
  n <- length(obs)

jumptimes <- c(0.00, 2.00, 2.50, 2.60, 4.80, 5.00, 8.00, 9.00, 13.50, 17.00, 20.48)
rates <- c(0.28308077, 0.35567019, 0.44679211, 0.52707776, 0.85724398, 0.64356224, 0.66275347, 0.49205892, 0.82856927, 0.54490179, 0.04882812)


 pointer <- c(1.8,5.3,7.9,9.3,15,17)
 piecerate <- c(1,10,2,8,2,10,2)
 U <- Uniform(jumptimes, rates, starting, ending, pointer, piecerate)

 index <- which(diff(U$U) <= 0.01)
 index
 Uniformtimes <- U$U
 Uniformtimes <- Uniformtimes[-index]
 JumpProb <- U$JumpProb[-index,]
 

 alpha <- 1.1
 zeta <- 2
 Var <- 1
 muprior <- 0
 sigmaprior <- 1
 Sigma <- 1
 

theta <- c(-0.526, 0.805, -0.852, 0.221, -0.579, 0.405, -0.571, -0.086, 1.096, 0.678, -0.405)

################################################################################
################################################################################
################################################################################
  CP <- res$CP
  Q <- res$Q
  theta <- res$theta
  changepoints <- Changepoints <- intensity <- jumprates <- CPnum <- NULL
################################################################################
################################################################################
###########################    Conjugate Priors   ##############################
################################################################################
################################################################################
  for (j in 1:iter){
  if (res$Number != 1){  
  rates <- -diag(Q)
  rates <- rates[-length(rates)]
  rates <- c(rates, 1/(ending-starting))
  flag <- FALSE
  while (flag==FALSE){
  A <- Uniform(CP, rates, starting, ending, pointer, piecerate)
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
  A <- Uniform(c(starting, ending), 1/(ending-starting), starting, ending, pointer, piecerate)
  Uniformtimes <- A$U
  index <- which(diff(c(starting,Uniformtimes,ending)) <= diff(eventtimes[1:2]))
  Uniformtimes <- Uniformtimes[-index]
  JumpProb <- A$JumpProb[-index, ]
  if ((length(Uniformtimes) >=3) && (length(Uniformtimes) <= 300)) flag <- TRUE
  }
  }
  print(j)
  
  res <- FFBS.Norm(Uniformtimes, JumpProb, eventtimes, obs, theta)
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
  if (j > burnin){
  changepoints <- c(changepoints, list(cp))
  Changepoints <- c(Changepoints, list(CP))
  jumprates <- c(jumprates, list(Q))
  CPnum <- c(CPnum, num)
  intensity <- c(intensity, list(theta))
        }
  }      
################################################################################
################################################################################
################################################################################
  CP <- res$CP
  Q <- res$Q
  theta <- res$theta
  Var <- res$Var
  changepoints <- Changepoints <- intensity <- Variance <- jumprates <- CPnum <- NULL
################################################################################
################################################################################
###########################Objective Priors##############################
################################################################################
################################################################################

  for (j in 1:iter){
  if (res$Number != 1){  
  rates <- -diag(Q)
  rates <- rates[-length(rates)]
  rates <- c(rates, 1/(ending-starting))
  flag <- FALSE
  while (flag==FALSE){
  A <- Uniform(CP, rates, starting, ending, pointer, piecerate)
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
  A <- Uniform(c(starting, ending), 1/(ending-starting), starting, ending, pointer, piecerate)
  Uniformtimes <- A$U
  index <- which(diff(c(starting,Uniformtimes,ending)) <= diff(eventtimes[1:2]))
  Uniformtimes <- Uniformtimes[-index]
  JumpProb <- A$JumpProb[-index, ]
  if ((length(Uniformtimes) >=3) && (length(Uniformtimes) <= 300)) flag <- TRUE
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
  
################################################################################
################################################################################
################################################################################
  CP <- res$CP
  Q <- res$Q
  theta <- res$theta
  Var <- res$Var
  changepoints <- Changepoints <- intensity <- Variance <- jumprates <- CPnum <- NULL
################################################################################
################################################################################
############################  Monte Carlo EM  ##################################
################################################################################
################################################################################

  for (j in 1:iter){
  if (res$Number != 1){  
  rates <- -diag(Q)
  rates <- rates[-length(rates)]
  rates <- c(rates, 1/(ending-starting))
  flag <- FALSE
  while (flag==FALSE){
  A <- Uniform(CP, rates, starting, ending, pointer, piecerate)
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
  A <- Uniform(c(starting, ending), 1/(ending-starting), starting, ending, pointer, piecerate)
  Uniformtimes <- A$U
  index <- which(diff(c(starting,Uniformtimes,ending)) <= diff(eventtimes[1:2]))
  Uniformtimes <- Uniformtimes[-index]
  JumpProb <- A$JumpProb[-index, ]
  if ((length(Uniformtimes) >=3) && (length(Uniformtimes) <= 300)) flag <- TRUE
  }
  }
  print(j)
 
  res <- MCEM.Objective(Uniformtimes, JumpProb, eventtimes, obs, theta, Var)
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
             