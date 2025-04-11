library(LaplacesDemon)

FFBS.Objective <- function(Uniformtimes, JumpProb, eventtimes, obs, param, Var){
n <- length(obs)
k <- length(Uniformtimes)
filtering <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
filtering[1,] <- c(1, rep(0,k))

JumpProb <- rbind(JumpProb, JumpProb[k,])
Uniformtimes <- c(Uniformtimes, ending)

for (i in 2:(k+1)){
   for (j in 1:(i-1))  filtering[i,j+1] = filtering[i-1, j+1]*JumpProb[i,1] + filtering[i-1, j]*JumpProb[i,2]
   filtering[i,1] = filtering[i-1,1]*JumpProb[i,1]
   ff <- logLik.Norm(eventtimes, obs, param, Var, Uniformtimes[i-1], Uniformtimes[i])
   filtering[i,] = log(filtering[i,]) + ff
   filtering[i,] = exp(filtering[i,] - max(filtering[i,]))
   filtering[i,] <- filtering[i,]/sum(filtering[i,])
   }
######  Sampling latent states  ###### 
S <- rep(1, k+1)
S[k+1] = min(sample(k+1, size=2, prob=filtering[k+1,], replace=T))
if(S[k+1] == 1) S[1:k] = 1 else{
for (i in k:2){
  Prob <- rep(0, 2)
  active = c(S[i+1]-1, S[i+1])
  Prob[2] = filtering[i, active[2]] * JumpProb[i+1, 1]
  Prob[1] = filtering[i, active[1]] * JumpProb[i+1, 2]
  S[i] = min(sample(active, size=2, prob=Prob, replace=T))
  if(S[i]==1){
              S[1:(i-1)] = 1
              break
             }
  }
  }
######  Extract CP     #######
CP <- NULL
CP <- Uniformtimes[which(diff(S)==1)]
m <- length(CP)
cp <- NULL
temp <- NULL
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
theta <- rep(NA, m+1)
#precision <- rep(NA, m+1)
for(i in 1:m){
Q[i, i] <- -rgamma(1, alpha + 1, zeta + sojourn[i])
Q[i, (i+1)] <- -Q[i,i]
t <- mean(obs[cp[i]:(cp[i+1]-1)])
theta[i] = rnorm(1, mean=t, sd=Var)
temp <- c(temp, obs[cp[i]:(cp[i+1]-1)]-theta[i])
}
t <- mean(obs[cp[m]:n])
theta[m+1] = rnorm(1, mean=t, sd=Var)
temp <- c(temp, obs[cp[m]:n]-theta[m+1])
cp <- cp[-1]
} else{
       cp <- NULL
       Q <- NULL
       t <- mean(obs)
       theta <- rnorm(1, mean=t, sd=Var)
       temp <- obs-theta
      }
CP <- c(0, CP, eventtimes[n])
Var <- sum(temp^2)/(n-1)
Var <- sqrt(rinvchisq(1, df=(n-1), scale=Var))
cat("\n Number:", m, "\n")
cat("Locations:", cp, "\n")
cat("Parameters:", theta, "\n")  
return(list(Number=m+1, cp=cp, CP=CP, Q=Q, states= S, theta=theta, Var=Var, filtering=filtering))     
}