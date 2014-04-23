library(rjags)
set.seed(201404231)
#############modelm - using point mass mixture prior for signals###############
modelm <- "
model{
# likelihood 
for (i in 1:length(y)){
y[i] ~ dpois(lambda[i])
log(lambda[i]) <- alpha[gene[i]] + (-1)^line[i]*tau[gene[i]] + beta[gene[i]]*cov[i]   
}
# prior level 1
for (i in 1:ngene){
alpha[i] ~ dnorm(0,1/10000)

tau[i] <- (1-bintau[i])*normtau[i]
bintau[i] ~ dbern(pitau)
normtau[i] ~ dnorm(0,1/sigmatau^2)
beta[i] <- (1-binbeta[i])*normbeta[i]
binbeta[i] ~ dbern(pibeta)
normbeta[i] ~ dnorm(0,1/sigmabeta^2)
}
#prior level 2
pitau ~ dbeta(8,1)
pibeta ~ dbeta(8,1)
sigmatau ~ dunif(0,100)
sigmabeta ~ dunif(0,100)
}
"

#########modelh - using horseshoe prior for the signals##############
modelh <- "
model{
# likelihood 
for (i in 1:length(y)){
  y[i] ~ dpois(lambda[i])
  log(lambda[i]) <- alpha[gene[i]] + (-1)^line[i]*tau[gene[i]] + beta[gene[i]]*cov[i]   
}
# prior level 1
for (j in 1:ngene){
  alpha[j] ~ dnorm(0,1/10000)
  tau[j] ~ dnorm(0, 1/sigmatauj[j]^2)
  sigmatauj[j] ~ dt(0, 1/sigmatau^2, 1) T(0,)
  beta[j] ~ dnorm(0, 1/sigmabetaj[j]^2)
  sigmabetaj[j] ~ dt(0, 1/sigmabeta^2, 1) T(0,)
}
# prior level 2
sigmatau ~ dt(0, 1, 1) T(0,)
sigmabeta ~ dt(0, 1, 1) T(0,)
}
"

####Simulation data ################
library(reshape)
K <- 12
x <- rnorm(2*K, 0,1)
ngene <- 100
# prior level 2
mualpha <- 3
sigmaalpha <- 2

pitau <-  0.8
mutau <- .5
sigmatau <- 0.25

pibeta <- 0.8
mubeta <- .5
sigmabeta  <- 0.25
# prior level 1
bintau <- rbinom(ngene,1,pitau)
tau <-(1-bintau)*rnorm(ngene, mutau, sigmatau)
beta <- (1-rbinom(ngene,1,pibeta))*rnorm(ngene, mubeta, sigmabeta)
alpha <- rnorm(ngene, mualpha,sigmaalpha )
lambda <- matrix(0, ncol = 2*K, nrow = ngene)
count <- matrix(0, ncol = 2*K, nrow = ngene)
for (j in 1:ngene){
  for (k in 1:K){
    lambda[j,k] <- exp(alpha[j] - tau[j] + beta[j]*x[k])
    lambda[j, k+K] <- exp(alpha[j] + tau[j] + beta[j]*x[k+K])
    count[j,k] <- rpois(1, lambda[j,k])
    count[j,k+K] <- rpois(1, lambda[j,k+K])
  }
}
melt_count <- melt(count)
melt_count$line <- NULL
melt_count$line[melt_count$X2%in%c(1:K)] <- 1
melt_count$line[melt_count$X2%in%c((K+1):(2*K))] <- 2
melt_count$cov <- NULL
for(i in 1:(2*K)) melt_count$cov[melt_count$X2==i] <- x[i]
data <- list(y = melt_count$value, 
             gene = melt_count$X1, 
             line = melt_count$line, 
             cov = melt_count$cov, 
             ngene = ngene)


######### run_mm_simulationdata #################
m0 <- proc.time()
mm <- jags.model(textConnection(modelm), data,n.chains = 1) # mix point mass
resm <- coda.samples(mm, c("tau","alpha","beta","pitau","pibeta", 
                           "binbeta",
                           "bintau", "sigmatau","sigmabeta"), 5000) # mix point mass
proc.time() -m0

######## plot_signal_tau_mm_simdata ############
which(tau!=0)
tau[which(tau!=0)]
mean(tau!=0)
plot(resm[, paste("tau[", which(tau!=0),"]",sep ="")])
plot(resm[, paste("bintau[", which(tau!=0),"]",sep ="")])
tau_est <- which(apply(resm[[1]][,paste("bintau[", 1:ngene,"]",sep ="")], 2, 
      function(x) mean(abs(1-x))) > 0.5)
tau_signal <- which(tau!=0)
sum(tau_est %in% tau_signal)
length(tau_signal)
length(tau_est) - sum(tau_est %in% tau_signal)

length(tau_est)
length(tau_signal)
apply(resm[[1]][,paste("binbeta[", 1:ngene,"]",sep ="")], 2, 
      function(x) mean(abs(1-x)))
mean(beta!=0)
which(beta!=0)
tau[which(beta!=0)]

which(apply(resm[[1]][,paste("beta[", 1:ngene,"]",sep ="")], 2, 
      function(x) mean(abs(x) >0))>0.5)

######## nonsignal_tau_mm_simdata #############
which(tau==0)
plot(resm[, paste("tau[", which(tau==0),"]",sep ="")])

apply(resm[[1]][,paste("tau[", 1:ngene,"]",sep ="")], 2, 
      function(x) mean(abs(x)>.1))


######## plot_signal_beta_mm_simdata ##########
which(beta!=0)
beta[which(beta!=0)]
plot(resm[, paste("beta[", which(beta!=0),"]",sep ="")])


######## nonsignal_beta_mm_simdata ############
which(beta==0)
plot(resm[, paste("beta[", which(beta==0),"]",sep ="")])

apply(resm[[1]][,paste("beta[", 1:ngene,"]",sep ="")], 2, 
      function(x) mean(abs(x)>0.1))
beta
####### plot alpha_mm_simdata ################
# alpha
# plot(resm[, paste("alpha[", 1:ngene,"]",sep ="")])
# 
# apply(resm[[1]][,paste("alpha[", 1:ngene,"]",sep ="")], 2, 
#       function(x) mean(abs(x)>1.5))
# 

####### plot sigmatau_mm_simdata ################

plot(resm[,"pibeta"])
beta
######### run_mh_simulationdata  #############

m0 <- proc.time()
mh <- jags.model(textConnection(modelh), data,n.chains = 1)
resh <- coda.samples(mh, c("tau","alpha","beta",
                          "sigmatauj", "sigmabetaj",  "sigmatau","sigmabeta"), 10000) 
proc.time() -m0

apply(resh[[1]][,paste("sigmatauj[", 1:ngene,"]",sep ="")], 2, 
      function(x) mean(1-1/(1+x^2)))

apply(resh[[1]][,paste("sigmabetaj[", 1:ngene,"]",sep ="")], 2, 
      function(x) mean(1-1/(1+x^2)))

######### plot_signal_tau_mh_simdata ##########
which(tau!=0)
plot(resh[, paste("tau[", which(tau!=0),"]",sep ="")])
tau[which(tau!=0)]

acf(resh[[1]][,"tau[6]"])
######## nonsignal_tau_mh_simdata ############
which(tau==0)
plot(resh[, paste("tau[", which(tau==0),"]",sep ="")])

apply(resh[[1]][,paste("tau[", 1:ngene,"]",sep ="")], 2, 
      function(x) mean(abs(x)>.1))

########## plot_signal_beta_mh_simdata #######

which(beta!=0)
beta[which(beta!=0)]
plot(resh[, paste("beta[", which(beta!=0),"]",sep ="")])


########## nonsignal_beta_mh_simdata #########
which(beta==0)
plot(resh[, paste("beta[", which(beta==0),"]",sep ="")])

apply(resh[[1]][,paste("beta[", 1:ngene,"]",sep ="")], 2, 
      function(x) mean(abs(x)>.1))

####### plot alpha_mh_simdata ################
alpha
plot(resh[, paste("alpha[", 1:ngene,"]",sep ="")])

apply(resh[[1]][,paste("alpha[", 1:ngene,"]",sep ="")], 2, 
      function(x) mean(abs(x)>.5))

