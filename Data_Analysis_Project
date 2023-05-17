library(tidyverse)
library(R2jags)
# Set directory
setwd("~/Documents/234")
# Load the data
load("dap1data.RData")

# Add predictor indicators
nursebp_final <- nursebp %>%
  mutate(famhist2 = ifelse(famhist == 2, 1, 0)) %>%
  mutate(famhist = ifelse(famhist == 1, 1, 0))  %>%
  mutate(famhistxwork = famhist*work) %>%
  mutate(famhist2xwork = famhist2*work)
nursebp_final <- nursebp_final[,c(1, 2, 3, 6, 4, 5, 7, 8)]

# Make count vector for model
nursebp_count <- nursebp_final %>%
  group_by(idnum) %>%
  count()
summary(nursebp_count)
nursebp_count_vector <- as.vector(t(nursebp_count [,2]))

# Make y matrix for model
test <- nursebp_final %>%
  group_by(idnum) %>%
  summarise(SYS = list(SYS)) %>%
  select(-idnum) %>%
  t()
max.leng <- max(sapply(test, length))
corrected_list <- lapply(test, function(x) {c(x,rep(NA, max.leng-length(x)))})
mat <- do.call(rbind, corrected_list)
final <- as.matrix(mat)

# Make cumulative count vector for model
cumcount <- c(rep(0,203))
for (i in 1:203){
  cumcount[i] <- sum(head(nursebp_count_vector, i-1))
}

# Summary
mysummary = function(invector) {
  c(mean(invector), sd(invector), quantile(invector, .025), 
    quantile(invector,.975),
    length(invector[invector<0])/length(invector),
    length(invector[invector==0])/length(invector),
    length(invector[invector>0])/length(invector))
}
# Create the model
sink("bloodpressure.txt")
cat("
model
{
  for(i in 1:203=
  {
    for(j in 1:obscount[i]){
      s[i,j] <- cumcount[i] + j
  	  y[i,j] ~ dnorm( mu[i,j], tau.epsilon) 
  	  mu[i,j] <- inprod(x[s[i,j],],beta[]) + gamma[i]
    }
    gamma[i] ~ dnorm(0, tau.gamma)
  }
  for(k in 1:7) {
  	beta[k] ~ dnorm( m[k], prec[k] )
  }
  tau.gamma ~ dgamma(da,db) 
  tau.epsilon ~ dgamma(ta,tb)
  sigma <- 1/sqrt(tau.epsilon)
  sqrtd <- 1/sqrt(tau.gamma)
}
", fill = TRUE)
sink()

# Declare all vectors, matrices, and priors
bpdata = list(y = final, cumcount = cumcount,
  x = matrix(
    data=c(rep(1,9573), nursebp_final[,3], nursebp_final[,4], nursebp_final[,5], 
      nursebp_final[,6], nursebp_final[,7], nursebp_final[,8]), 
    byrow=F, ncol=7), 
  obscount = nursebp_count_vector, da=61, db = 5184, ta=237.25, tb= 14066.525, 
  m=c(114.0,3.2, 5.3, 2.6, 0.22, 0, 0), prec = c(
    0.001,0.20,0.1, 1, 1.1, 0.025, 0.025))

# Initial values and parameters
bpinits = rep(list(list( beta = c(117.0,0,0,0,0,0,0), 
  gamma=as.vector(rnorm(203, mean = 0, sd = 7)), tau.gamma=1, 
  tau.epsilon = 1)), 5)

bpparameters = c("beta","tau.epsilon", "tau.gamma", "sqrtd", "sigma")

# Run the model
run1 = jags(bpdata, bpinits, bpparameters,
               "bloodpressure.txt", n.chains=5, n.iter=10000, 
               n.burnin=1000, n.thin=1)
print(run1)
# Table 1
temp <- run1$BUGSoutput$sims.matrix
out <- t(apply(temp[,1:11], 2, mysummary))
colnames(out) <- c("mean", "sd", " 2.50%", "97.5%", "P(<0|Y)", "P(=0|Y)",
  "P(>0|Y)")
out
out_final <- round(out,2)
write.table(out_final, file = "Out.txt", sep = ",", quote = FALSE, 
  row.names = F)
# Posterior Plots
xaxisvals = seq(from = -10,to = 15, by = 0.01)
## Figure 1
par(mfrow=c(1,2))
par(mar=c(3.1,4.1,2.1,2.1))
plot(density(temp[,9]), main="sigma", xlim=c(12, 15))
abline(v=0)
plot(density(temp[,10]), main="sqrtD", xlim=c(4,12))
abline(v=0)
xaxisvals = seq(from = -10,to = 15, by = 0.01)
xaxisvals1 = seq(from = 70,to = 180, by = 0.01)
## Figure 2
par(mfrow=c(3,3))
par(mar=c(3.1,4.1,2.1,2.1))
plot(density(temp[,1]), main="beta[0]", xlim = c(110, 130), col = "red")
lines(xaxisvals1, dnorm(xaxisvals1, mean = 114, sd = sqrt(1000)), col = "blue")
plot(density(temp[,2]), main="beta[1]", col = "red", xlim = c (-5,10))
lines(xaxisvals, dnorm(xaxisvals, mean = 3.2, sd = 2.23), col = "blue")
abline(v=0)
plot(density(temp[,3]), main="beta[2]", col = "red", xlim = c (-5,15))
lines(xaxisvals, dnorm(xaxisvals, mean = 5.3, sd = sqrt(10)), col = "blue")
abline(v=0)
plot(density(temp[,4]), main="beta[3]",col = "red", xlim = c (-2,7))
lines(xaxisvals, dnorm(xaxisvals, mean = 2.6, sd = 1), col = "blue")
abline(v=0)
plot(density(temp[,5]), main="beta[4]", col = "red",)
lines(xaxisvals, dnorm(xaxisvals, mean = 0.22, sd = sqrt(1/1.1)), col = "blue")
abline(v=0)
plot(density(temp[,6]), main="beta[5]", col = "red", ylim = c(0, 1.0))
lines(xaxisvals, dnorm(xaxisvals, mean = 0, sd = sqrt(40)), col = "blue")
abline(v=0)
plot(density(temp[,7]), main="beta[6]", col = "red", ylim = c(0, 1.0))
lines(xaxisvals, dnorm(xaxisvals, mean = 0, sd = sqrt(40)), col = "blue")
abline(v=0)
legend(x = "right", legend=c("Prior", "Posterior"), col=c("blue", "red"), lty=1,
  lwd=1, cex = 0.8, seg.len = 1, inset = c(-3,0), xpd = NA, 
  title = "Distribution")

# Convergence
## Autocorrelation plots
temp2= run1$BUGSoutput$sims.array
dim(temp2)
par(mfrow=c(3,3))     
par(mar=c(3.1,4.1,2.1,2.1))   
acf(temp2[,1,1], main="")     
mtext("beta[0]",side=3, line=1, cex=.8)
acf(temp2[,1,2], main="") 
mtext("beta[1]",side=3, line=1, cex=.8)
acf(temp2[,1,3], main="") 
mtext("beta[2]",side=3, line=1, cex=.8)
acf(temp2[,1,4], main="") 
mtext("beta[3]",side=3, line=1, cex=.8)
acf(temp2[,1,5], main="") 
mtext("beta[4]",side=3, line=1, cex=.8)
acf(temp2[,1,6], main="") 
mtext("beta[5]",side=3, line=1, cex=.8)
acf(temp2[,1,7], main="") 
mtext("beta[6]",side=3, line=1, cex=.8)
acf(temp2[,1,8], main="") 
mtext("SqrtD",side=3, line=1, cex=.8)
acf(temp2[,1,10], main="") 
mtext("Sigma",side=3, line=1, cex=.8)
par(mfrow=c(3,3))
par(mar=c(4.1,4.1,2.1,2.1))

## Time-series plots
plot(1:900,temp2[1:900,1,1], type="l", xlab="iteration", ylab="")  
mtext("beta[0]", side=3, line=1, cex=.8)
plot(1:900,temp2[1:900,1,2], type="l", xlab="iteration", ylab="")  
mtext("beta[1]", side=3, line=1, cex=.8)
plot(1:900,temp2[1:900,1,3], type="l", xlab="iteration", ylab="")  
mtext("beta[2]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,4]),temp2[,1,4], type="l", xlab="iteration", ylab="")  
mtext("beta[3]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,5]),temp2[,1,4], type="l", xlab="iteration", ylab="")  
mtext("beta[4]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,6]),temp2[,1,4], type="l", xlab="iteration", ylab="")  
mtext("beta[5]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,7]),temp2[,1,4], type="l", xlab="iteration", ylab="")  
mtext("beta[6]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,8]),temp2[,1,4], type="l", xlab="iteration", ylab="")  
mtext("SqrtD", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,10]),temp2[,1,4], type="l", xlab="iteration", ylab="")  
mtext("Sigma", side=3, line=1, cex=.8)

#Sensitivity Analysis (Table 2)
## Low variance
bpdata2 = list(y = final, cumcount = cumcount,
  x = matrix(
    data=c(rep(1,9573), nursebp_final[,3], nursebp_final[,4], nursebp_final[,5], 
      nursebp_final[,6], nursebp_final[,7], nursebp_final[,8]), 
    byrow=F, ncol=7), 
  obscount = nursebp_count_vector, da=61, db = 5184, ta=237.25, tb= 14066.525, 
  m=c(114.0,3.2, 5.3, 2.6, 0.22, 0, 0), prec = c(
    0.003, 0.60, 0.3, 3, 3.3, 0.075, 0.075))
run2 = jags(bpdata2, bpinits, bpparameters,
  "bloodpressure.txt", n.chains=5, n.iter=5000, 
  n.burnin=1000, n.thin=1)
print(run2)
temp2 <- run2$BUGSoutput$sims.matrix
out2 <- t(apply(temp2[,1:11], 2, mysummary))
colnames(out2) <- c("mean", "sd", " 2.50%", "97.5%", "P(<0|Y)", "P(=0|Y)","P(>0|Y)")
out2_final <- round(out2,2)
write.table(out2_final, file = "Out2.txt", sep = ",", quote = FALSE, row.names = F)

## High Variance
bpdata3 = list(y = final, cumcount = cumcount,
  x = matrix(
    data=c(rep(1,9573), nursebp_final[,3], nursebp_final[,4], nursebp_final[,5], 
      nursebp_final[,6], nursebp_final[,7], nursebp_final[,8]), 
    byrow=F, ncol=7), 
  obscount = nursebp_count_vector, da=61, db = 5184, ta=237.25, tb= 14066.525,
  m=c(114.0,3.2, 5.3, 2.6, 0.22, 0, 0), prec = c(
    0.001/3,0.20/3, 0.1/3, 1/3, 1.1/3, 0.025/3, 0.025/3))
run3 = jags(bpdata3, bpinits, bpparameters,
  "bloodpressure.txt", n.chains=5, n.iter=5000, 
  n.burnin=1000, n.thin=1)
print(run3)
temp3 <- run3$BUGSoutput$sims.matrix
out3 <- t(apply(temp3[,1:11], 2, mysummary))
colnames(out3) <- c("mean", "sd", " 2.50%", "97.5%", "P(<0|Y)", "P(=0|Y)","P(>0|Y)")
out3_final <- round(out3,2)
write.table(out3_final, file = "Out3.txt", sep = ",", quote = FALSE, row.names = F)
