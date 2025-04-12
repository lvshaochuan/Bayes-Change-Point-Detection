jumptimes <- c(0.00,10.50,15.00,16.50,18.20,20.00,23.50,24.30,25.00,25.50,27.50,39.57)
rates <- c(0.09185865,0.34486715,0.27146221,0.51170168,0.54673125,0.35096016,0.33,0.9,1.36,0.6,0.025)
U <- Uniform(jumptimes, rates)
Uniformtimes <- U$U
JumpProb <- U$JumpProb
index <- which(diff(U$U)<=0.01)
Uniformtimes <- Uniformtimes[-index]
JumpProb <- JumpProb[-index,]
length(obs)
[1] 3957
eventtimes <- seq(0,39.57,by=0.01)
eventtimes <- eventtimes[-1]
obs <- (obs-115000)/10000
mean(obs)
[1] 0.1760483
sd(obs)
[1] 0.8109748
param <- c( -0.5,  1.5 , 2.0,  0.0 , 1.5,  0.5,  2.0 , 0.5  ,1.5 , 0.0,-0.5)
> alpha
[1] 1.1
> zeta
[1] 2
> Var <- 1
> res <- MCEM.Objective(Uniformtimes, JumpProb, eventtimes, obs, param, Var)




################################# Well-log by FFBS.Norm2.R #####################
 position <- NULL
 for (i in 1:length(CPnum)){
 temp <- length(Changepoints[[i]])
 position <- c(position, Changepoints[[i]][-c(1,temp)])
 }
 

savpar <- par(mfrow=c(3,1))
plot(obs, xlab="", ylab="", xlim=c(0, 4050), main="(a) Well-log Data", cex=0.7)
#abline(v=c(1046,1486,1645,1825,2006,2368,2428,2490,2550,2727), lty="dashed", col="red")
hist(position/eventtimes[1], 800, xlim=c(0, 4050), xlab="Changepoint Location", main="(b) Changepoint Locations")
box()
hist(CPnum[1001:length(CPnum)]-1, 80, xlab="Number of Changepoints", main="(c) Number of Changepoints", col="grey")
#hist(parameter, 500, xlab=expression(mu), main="(c) Histogram")
box()
par(savpar) 

hist(CPnum[1:length(CPnum)]-1, 800, xlab="Number of Changepoints", main="(c) Number of Changepoints", col="grey")

savpar <- par(mfrow=c(2,1))
plot(obs, xlab="", ylab="", xlim=c(0, 3957), main="(a) Well-log Data", cex=0.5)
abline(v=c(1046,1486,1645,1825,2006,2368,2428,2490,2550,2727), lty="dashed", col="red")
hist(position/eventtimes[1], 1500, xlim=c(0, 3957), xlab="Changepoint Location", main="(b) Changepoint Locations")
box()
par(savpar)


temp <- NULL
for (i in 1176:2175) temp <- rbind(temp, intensity[[i]])
postmean <- apply(temp, 2, mean)

Qrate <- jumprates[[1176]]
for (i in 1177:2175){
Qrate <- Qrate + jumprates[[i]]
}
Qrate <- Qrate/1000

path <- Viterbi(diff(eventtimes), obs[-1], Qrate, delta=c(1,rep(0,10)), postmean, rep(1,11))
which(diff(path)==1)
 [1] 1046 1486 1645 1825 2006 2368 2428 2490 2550 2727


 ########################### Well-log by FFBS.new.Objective.R ##################################
 parameter <- NULL
 for (i in 1:length(CPnum)){
 temp <- length(Changepoints[[i]])
 parameter <- c(parameter, intensity[[i]][-c(1,temp)])
 } 
 
split.screen(c(2,1))
split.screen(c(1,2),screen=1)

screen(3)
hist(CPnum[1001:3000]-1, 500, xlab="Number", main="(a) Number of Changepoints", col="grey")
box()                                                                                                 
screen(4)
#hist(thetas, 300, xlab=expression(mu), main="(b) Histogram", col="grey")
hist(parameter, 300, xlab=expression(mu), main="(b) Histogram", col="grey")
box()

screen(2)
#hist(position, 200, xlab="Changepoint Location", main=" (c)Changepoint Locations")
hist(position*100, 1000, xlab="Changepoint Location", main="(c) Changepoint Locations")
box()

############################## Well-log by MCEMnewRun.Objective.R##################################################

split.screen(c(2,1))
split.screen(c(1,2),screen=1)

screen(3)
hist(parameter, 500, xlab=expression(mu), main="(a) Histogram", col="grey")
box()                                                                                                 
screen(4)
hist(Variance, 300, xlab=expression(sigma^2), main="(b) Histogram", col="grey")
box()

screen(2)
#hist(position, 200, xlab="Changepoint Location", main=" (c)Changepoint Locations")
hist(position*100, 1000, xlab="Changepoint Location", main="(c) Changepoint Locations")
box()

 x11()
 split.screen(c(2,1))
 split.screen(c(1,2),screen=1)
 screen(3)
 hist(parameter, 500, xlab=expression(mu), main="(a) Histogram", col="grey")
 box()                                                                                                 
 screen(4)
 hist(sqrt(Variance), 300, xlab=expression(sigma), main="(b) Histogram", col="grey")
 box()
 
 screen(2)
 #hist(position, 200, xlab="Changepoint Location", main=" (c)Changepoint Locations")
 hist(position*100, 1000, xlab="Changepoint Location", main="(c) Changepoint Locations")
 box()

############################## Well-log by Depend.Run.R##################################################

split.screen(c(2,1))
split.screen(c(1,2),screen=1)

screen(3)
plot(obs, xlab="", ylab="", xlim=c(0, 3957), main="(a) Well-log Data", cex=0.5)
#abline(v=c(1046,1486,1645,1825,2006,2368,2428,2490,2550,2727), lty="dashed", col="red")                                                                                       
screen(4)
#hist(thetas, 300, xlab=expression(mu), main="(b) Histogram", col="grey")
hist(parameter, 300, xlab=expression(mu), main="(b) Histogram", col="grey")
box()

screen(2)
#hist(position, 200, xlab="Changepoint Location", main=" (c)Changepoint Locations")
hist(position*100, 1000, xlab="Changepoint Location", main="(c) Changepoint Locations")
box()


savpar <- par(mfrow=c(2,1))
#hist(position, 200, xlab="Changepoint Location", main=" (c)Changepoint Locations")
hist(position*100, 3000, xlab="Changepoint Location", main="(a) Changepoint Locations")
box()

#hist(thetas, 300, xlab=expression(mu), main="(b) Histogram", col="grey")
hist(parameter, 300, xlab=expression(mu), main="(b) Histogram ", col="grey")
box()
par(savpar)
