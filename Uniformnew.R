starting
ending
pt1  <- 2.5
pt2  <- 10.5
rate1  <- 3
rate2  <- 2
rate3  <- 3

Q <- res$Q
rates <- -diag(Q)
rates <- rates[-length(rates)]
rates <- c(rates, 1/(ending-starting))

Uniform <- function(jumptimes, rates){
##### Simulate Uniformisation times #####
initial <- starting
Poissontimes <- NULL
tag <- NULL
while(initial <= ending){
if (initial < pt1) {
Poissontime <-  rexp(1, rate1)
initial <- initial + Poissontime
if(initial <= pt1) {
   Poissontimes <- c(Poissontimes, initial)
   tag <- c(tag, 1)
   } else initial <- pt1
} else{
    if (initial < pt2){
        Poissontime <-  rexp(1, rate2)
        initial <- initial + Poissontime        
        if(initial <= pt2) {Poissontimes <- c(Poissontimes, initial)
        tag <- c(tag, 2)
        } else initial <- pt2
    } else{
    Poissontime <-  rexp(1, rate3)
    initial <- initial + Poissontime
    if(initial <= ending) {Poissontimes <- c(Poissontimes, initial)
       tag <- c(tag, 3)} else break
         }
 }
}
##### Tagged Uniformisation times #####
i <- 2
j <- 1
k <- length(jumptimes)-2
if(k != 0){tag0 <- NULL
while(i < (k+2)){
   if(jumptimes[i] <= pt1) {tag0 <- c(tag0, 1)}  else{
   if(jumptimes[i] <= pt2) {tag0 <- c(tag0, 2)} else tag0 <- c(tag0, 3)
   }
   i <- i + 1
  }
} else tag0 <- NULL
Poissontimes <- c(Poissontimes, jumptimes[-c(1,k+2)])
tag <- c(tag, tag0)
temp <- order(Poissontimes)
tag <- tag[temp]
Poissontimes <- Poissontimes[temp]
n <- length(Poissontimes)
##### Transition Probability Matrix #####
m <- length(tag0)
delimit <- NULL
mark <- rep(1/n, n)
if(m != 0){
for(i in 1:m){
  j = which(Poissontimes==jumptimes[i+1])
  delimit <- c(delimit, j)
  }
mark[1:delimit[1]-1] = rates[1] 
if(m >= 2) {for(i in 2:m) mark[delimit[i-1]:(delimit[i]-1)] = rates[i]}
#mark[delimit[m]:n] = rates[m + 1]
}
Rate <- rep(NA, n)
Rate[tag==1] = rate1
Rate[tag==2] = rate2
Rate[tag==3] = rate3
JumpProb <- matrix(rep(NA, 2*n), ncol=2)
for(i in 1:n){
JumpProb[i,2] = mark[i]/(Rate[i] + mark[i])
JumpProb[i,1] = 1- JumpProb[i,2]
}
return(list(U=Poissontimes, JumpProb=JumpProb))
}