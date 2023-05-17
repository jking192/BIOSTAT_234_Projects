library(R2jags)
library(tidyverse)
# Set directory
setwd("~/Documents/234")
# Load the data
data <- haven::read_sas("organized_data.sas7bdat")
#Import Data
data <- data[,c(7, 6, 10, 11, 3, 2, 4, 5, 9, 8, 1)]
y <- data$obesity 
x <- as.matrix(data[,2:11]) 
# Model
sink("FDAP_model.txt")
cat("
model 
{
  for(i in 1:N)
  {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- beta0 + inprod(x[i,],beta[]) 
  }
  beta0 ~ dnorm(mbeta0, precbeta0)
  baseodds <- exp(beta0)
  for(j in 1:10)
  {
    beta[j] ~ dnorm(m[j], prec[j])
    OR[j] <- exp(beta[j])
  }
}
", fill = TRUE)
sink()
# Summary
mysummary = function(invector) {
  c(mean(invector), sd(invector), quantile(invector, .025), 
    quantile(invector,.975),
    length(invector[invector<1])/length(invector),
    length(invector[invector>1])/length(invector))
}
# Declare priors and data
fdapdata <- list(y = y,
  x = x,
  m=c(-0.253,-0.15, -0.31, 0.34, -0.051, 0.051, 0.13, 0.50, 0.50, -0.089), 
  prec = c( 4.13, 2.40, 1.82, 0.45, 0.70, 0.94, 1.20, 1.17, 0.89, 6.56), 
  mbeta0 = -1.68, precbeta0 = 0.318, N = 6818)
# Initial values and parameters to track
fdapinits <- rep(list(list(beta0 = 0, beta = c(0,0,0,0,0,0,0,0,0,0))), 5)
fdapparameters <- c("beta0", "beta","OR", "baseodds")
# Run the model
run1 <- jags(fdapdata, fdapinits, fdapparameters,
            "FDAP_model.txt", n.chains=5, n.iter=8000, 
            n.burnin=1000, n.thin=1)
print(run1)

temp <- run1$BUGSoutput$sims.matrix
out <- t(apply(temp[,c(1:11)], 2, mysummary))
colnames(out) <- c("mean", "sd", "2.50%", "97.5%", "P < 1", 
                   "P > 1")
out_final <- round(out,2)
write.table(out_final, file = "Out.txt", sep = ",", quote = FALSE, 
            row.names = F)
# Posterior Plots
xaxisvals = seq(from = -5,to = 5, by = 0.01)
## Figure 2
par(mfrow=c(4,3))
par(mar=c(3.1,4.1,1,1))
plot(density((temp[,22])), main="Intercept",  col = "red", xlim = c(-3,3))
lines(xaxisvals, dnorm(xaxisvals, mean = -1.68, sd = 1.77), col = "blue")
abline(v=0)
plot(density((temp[,12])), main="Female", col = "red", xlim = c (-1,1))
lines(xaxisvals, dnorm(xaxisvals, mean = -0.253, sd = sqrt(1/4.13)), col = "blue")
abline(v=0)
plot(density(temp[,13]), main="Middleage", col = "red", xlim = c(-1,1))
lines(xaxisvals, dnorm(xaxisvals, mean = -0.15, sd = sqrt(1/2.40)), col = "blue")
abline(v=0)
plot(density(temp[,14]), main="Senior",col = "red", xlim = c(-1,1))
lines(xaxisvals, dnorm(xaxisvals, mean = -0.31, sd = sqrt(1/1.82)), col = "blue")
abline(v=0)
plot(density(temp[,15]), main="Poor", col = "red", xlim = c(-1,1))
lines(xaxisvals, dnorm(xaxisvals, mean = 0.34, sd = sqrt(1/0.45)), col = "blue")
abline(v=0)
plot(density(temp[,16]), main="Rich", col = "red", xlim = c(-2,2))
lines(xaxisvals, dnorm(xaxisvals, mean = -0.051, sd = sqrt(1/0.70)), col = "blue")
abline(v=0)
plot(density(temp[,17]), main="Highorless", col = "red", xlim = c(-5,5))
lines(xaxisvals, dnorm(xaxisvals, mean = 0.051, sd = sqrt(1/0.94)), col = "blue")
abline(v=0)
plot(density(temp[,18]), main="Somecollege", col = "red", xlim = c(-1,1))
lines(xaxisvals, dnorm(xaxisvals, mean = 0.13, sd = sqrt(1/1.20)), col = "blue")
abline(v=0)
plot(density(temp[,19]), main="Married", col = "red", xlim = c(-1,1))
lines(xaxisvals, dnorm(xaxisvals, mean = 0.50, sd = sqrt(1/1.17)), col = "blue")
abline(v=0)
plot(density(temp[,20]), main="Othermarital", col = "red", xlim = c(-1,1))
lines(xaxisvals, dnorm(xaxisvals, mean = 0.50, sd = sqrt(1/0.89)), col = "blue")
abline(v=0)
plot(density(temp[,21]), main="Smoking", col = "red",xlim = c(-1,1))
lines(xaxisvals, dnorm(xaxisvals, mean = -0.089, sd = sqrt(1/6.56)), col = "blue")
abline(v=0)
legend(x = "right", legend=c("Prior", "Posterior"), col=c("blue", "red"), lty=1,
       lwd=1, cex = 0.8, seg.len = 1, inset = c(-1.7,0), xpd = NA, 
       title = "Distribution")

# Convergence
## Autocorrelation plots
temp2= run1$BUGSoutput$sims.array
dim(temp2)
par(mfrow=c(4,3))     
par(mar=c(3.1,4.1,2.1,2.1))   
acf(temp2[,1,22], main="")     
mtext("beta[0]",side=3, line=1, cex=.8)
acf(temp2[,1,12], main="") 
mtext("beta[1]",side=3, line=1, cex=.8)
acf(temp2[,1,13], main="") 
mtext("beta[2]",side=3, line=1, cex=.8)
acf(temp2[,1,14], main="") 
mtext("beta[3]",side=3, line=1, cex=.8)
acf(temp2[,1,15], main="") 
mtext("beta[4]",side=3, line=1, cex=.8)
acf(temp2[,1,16], main="") 
mtext("beta[5]",side=3, line=1, cex=.8)
acf(temp2[,1,17], main="") 
mtext("beta[6]",side=3, line=1, cex=.8)
acf(temp2[,1,18], main="") 
mtext("beta[7]",side=3, line=1, cex=.8)
acf(temp2[,1,19], main="") 
mtext("beta[8]",side=3, line=1, cex=.8)
acf(temp2[,1,20], main="") 
mtext("beta[9]",side=3, line=1, cex=.8)
acf(temp2[,1,21], main="") 
mtext("beta[10]",side=3, line=1, cex=.8)


## Time-series plots
par(mfrow=c(4,3))
par(mar=c(4.1,4.1,2.1,2.1))
plot(1:length(temp2[,1,22]),temp2[,1,22], type="l", xlab="iteration", ylab="")  
mtext("beta[0]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,12]),temp2[,1,12], type="l", xlab="iteration", ylab="")  
mtext("beta[1]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,13]),temp2[,1,13], type="l", xlab="iteration", ylab="")  
mtext("beta[2]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,14]),temp2[,1,14], type="l", xlab="iteration", ylab="")  
mtext("beta[3]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,15]),temp2[,1,15], type="l", xlab="iteration", ylab="")  
mtext("beta[4]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,16]),temp2[,1,16], type="l", xlab="iteration", ylab="")  
mtext("beta[5]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,17]),temp2[,1,17], type="l", xlab="iteration", ylab="")  
mtext("beta[6]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,18]),temp2[,1,18], type="l", xlab="iteration", ylab="")  
mtext("beta[7]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,19]),temp2[,1,19], type="l", xlab="iteration", ylab="")  
mtext("beta[8]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,20]),temp2[,1,20], type="l", xlab="iteration", ylab="")  
mtext("beta[9]", side=3, line=1, cex=.8)
plot(1:length(temp2[,1,21]),temp2[,1,21], type="l", xlab="iteration", ylab="")  
mtext("beta[10]", side=3, line=1, cex=.8)
# Sensitivity Analyses
## Very Low Variance
fdapdata2 <- list(y = y,
  x = x,
  m=c(-0.253,-0.15, -0.31, 0.34, -0.051, 0.051, 0.13, 0.50, 0.50, -0.089), 
  prec = c( 4.13*9, 2.40*9, 1.82*9, 0.45*9, 0.70*9, 
    0.94*9, 1.20*9, 1.17*9, 0.89*9, 6.56*9), 
  mbeta0 = -1.68, precbeta0 = 0.318, N = 6818)
run2 <- jags(fdapdata2, fdapinits, fdapparameters,
             "FDAP_model.txt", n.chains=5, n.iter=10000, 
             n.burnin=1000, n.thin=1)
## Low Variance
fdapdata3 <- list(y = y,
  x = x,
  m=c(-0.253,-0.15, -0.31, 0.34, -0.051, 0.051, 0.13, 0.50, 0.50, -0.089), 
  prec = c( 4.13*3, 2.40*3, 1.82*3, 0.45*3, 0.70*3, 0.94*3, 
    1.20*3, 1.17*3, 0.89*3, 6.56*3), 
mbeta0 = -1.68, precbeta0 = 0.318, N = 6818)
run3 <- jags(fdapdata3, fdapinits, fdapparameters,
             "FDAP_model.txt", n.chains=5, n.iter=10000, 
             n.burnin=1000, n.thin=1)
print(run2)
temp2 <- run2$BUGSoutput$sims.matrix
out2 <- t(apply(temp2[,c(1:11)], 2, mysummary))
colnames(out2) <- c("mean", "sd", "2.50%", "97.5%", "P < 1", 
                   "P > 1")
out2_final <- round(out2,2)
write.table(out2_final, file = "Out2_final.txt", sep = ",", quote = FALSE, 
            row.names = F)

print(run3)
temp3 <- run3$BUGSoutput$sims.matrix
out3 <- t(apply(temp3[,c(1:11)], 2, mysummary))
colnames(out3) <- c("mean", "sd", "2.50%", "97.5%", "P < 1", 
                    "P > 1")
out3_final <- round(out3,2)
write.table(out3_final, file = "Out3_final.txt", sep = ",", quote = FALSE, 
            row.names = F)

# T-distributed regression coefficients
# Model
sink("FDAP_model_t.txt")
cat("
model 
{
  for(i in 1:N)
  {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- beta0 + inprod(x[i,],beta[]) 
  }
  beta0 ~ dt(mbeta0, precbeta0, 5)
  baseodds <- exp(beta0)
  for(j in 1:10)
  {
    beta[j] ~ dt(m[j], prec[j], 5)
    OR[j] <- exp(beta[j])
  }
}
", fill = TRUE)
sink()

run4 <- jags(fdapdata, fdapinits, fdapparameters,
             "FDAP_model_t.txt", n.chains=5, n.iter=10000, 
             n.burnin=100, n.thin=1)
temp4 <- run4$BUGSoutput$sims.matrix
out4 <- t(apply(temp4[,c(1:11)], 2, mysummary))
colnames(out4) <- c("mean", "sd", "2.50%", "97.5%", "P < 1", 
                   "P > 1")
out4_final <- round(out4,2)
write.table(out4_final, file = "Out4_final.txt", sep = ",", quote = FALSE, 
            row.names = F)
