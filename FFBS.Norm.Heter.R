require(LaplacesDemon)
FFBS.Norm <- function(Uniformtimes, JumpProb, eventtimes, obs, param, Sigma){
n <- length(obs)
k <- length(Uniformtimes)
filtering <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
filtering[1,] <- c(1, rep(0,k))

JumpProb <- rbind(JumpProb, JumpProb[k,])
Uniformtimes <- c(Uniformtimes, ending)

for (i in 2:(k+1)){
   for (j in 1:(i-1))  filtering[i,j+1] = filtering[i-1, j+1]*JumpProb[i,1] + filtering[i-1, j]*JumpProb[i,2]
   filtering[i,1] = filtering[i-1,1]*JumpProb[i,1]
   ff <- logLik.Norm2(eventtimes, obs, param, Sigma, Uniformtimes[i-1], Uniformtimes[i])
   filtering[i,] = log(filtering[i,]) + ff
   filtering[i,] = exp(filtering[i,] - max(filtering[i,]))
   filtering[i,] <- filtering[i,]/sum(filtering[i,])
   }
######  Sampling latent states   
S <- rep(1, k+1)
S[k+1] = min(sample(k+1, size=5, prob=filtering[k+1,], replace=T))
if(S[k+1] == 1) S[1:k] = 1

for (i in k:2){
  Prob <- rep(0, 2)
  active = c(S[i+1]-1, S[i+1])
  scaling <- filtering[i, active[1]] + filtering[i, active[2]] 
  Prob[2] = filtering[i, active[2]]/scaling * JumpProb[i+1, 1]
  Prob[1] = filtering[i, active[1]]/scaling * JumpProb[i+1, 2]
  S[i] = min(sample(active, size=5, prob=Prob, replace=T))
  if(S[i]==1){
              S[1:(i-1)] = 1
              break
             }
  }
cat(S, "\n")
######  Extract CP 
CP <- NULL
CP <- Uniformtimes[which(diff(S)==1)]
m <- length(CP)
cp <- NULL
j <- 1
if(m > 0){
for (i in 1:n){
     if(eventtimes[i] >= CP[j]){
     cp <- c(cp, i)
     j <- j + 1
     if(j > length(CP)) break
     }
     }
N <- c(cp[1]-1, diff(cp), (n+1)-cp[m])
sojourn <- c(CP[1], diff(CP), eventtimes[n]-CP[m])
cp <- c(1, cp)
Q <- matrix(rep(0, (m+1)*(m+1)), nrow=m+1)
sigma2 <- rep(NA, m+1)
Mu <- rep(NA, m+1)

for(i in 1:m){
Q[i, i] <- -rgamma(1, alpha + 1, zeta + sojourn[i])
Q[i, (i+1)] <- -Q[i,i]
Xbar <- mean(obs[cp[i]:(cp[i+1]-1)])
KappaN <- Kappa0 + N[i]
MuN <- (Kappa0*Mu0 + sum(obs[cp[i]:(cp[i+1]-1)]))/KappaN
NvN <- Nv0 + N[i]
SigmaN <- (Nv0*Sigma0^2 + sum((obs[cp[i]:(cp[i+1]-1)] - Xbar)^2) + N[i]*Kappa0*(Mu0 - Xbar)^2/KappaN)/NvN  
sigma2[i] <- rinvchisq(1, df=NvN, scale=sqrt(SigmaN))
Mu[i] <- rnorm(1, mean=MuN, sd=sqrt(SigmaN/KappaN))      
}
Xbar <- mean(obs[cp[m]:n])
KappaN <- Kappa0 + N[m]
MuN <- (Kappa0*Mu0 + sum(obs[cp[m]:n]))/KappaN
NvN <- Nv0 + N[m]
SigmaN <- (Nv0*Sigma0^2 + sum((obs[cp[m]:n] - Xbar)^2) + N[m]*Kappa0*(Mu0 - Xbar)^2/KappaN)/NvN  
sigma2[m+1] <- rinvchisq(1, df=NvN, scale=sqrt(SigmaN))
Mu[m+1] <- rnorm(1, mean=MuN, sd=sqrt(SigmaN/KappaN)) 
#param[m] <- mean(obs[cp[m]:n])
#param[m] <- rnorm(1, mean(obs[cp[m]:n]), sd=sigma/N[m])
##Lambda[m] <- rgamma(1, a + N[m], b + sojourn[m])
cp <- cp[-1]
} else{
m <- 1
cp <- NULL
Q <- NULL
#param <- rnorm(1, mean(obs), sd=sigma/n)
Xbar <- mean(obs)
KappaN <- Kappa0 + n
MuN <- (Kappa0*Mu0 + sum(obs))/KappaN
NvN <- Nv0 + n
SigmaN <- (Nv0*Sigma0^2 + sum((obs - Xbar)^2) + n*Kappa0*(Mu0 - Xbar)^2/KappaN)/NvN  
sigma2 <- rinvchisq(1, df=NvN, scale=sqrt(SigmaN))
Mu <- rnorm(1, mean=MuN, sd=sqrt(SigmaN/KappaN))      

      }
CP <- c(0, CP, eventtimes[n])
cat("\n Number:", m, "\n")
cat("Locations:", cp, "\n")
cat("Parameters:", Mu, sigma2, "\n")  
return(list(Number=m, cp=cp, CP=CP, Q=Q, filtering=filtering, Sigma=sigma2, Mu=Mu))     
}