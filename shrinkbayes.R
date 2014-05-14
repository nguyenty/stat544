library(ShrinkBayes)
# load("data.RData")
# dat2 <- as.data.frame(cbind(data[[1]], data[[2]], data[[3]], data[[4]]))
# 
# 
# colnames(dat2) <- c("y", "gene", "line", "cov")
# 
# 
# library(reshape)
# dat3 <- cast(dat2, formula = gene ~ cov,  value= "y")
# str(dat3)
# dim(dat3)
# colnames(dat3)
# head(dat2)
# tail(dat2)
# count <- dat3[,-1]
# dim(count)
# head(count)
# cov <- as.numeric(colnames(count))
# line <- rep(c(1,2), each=12)
# line
# head(count)
# count[26,]
# count[81,]
# counts <- count
# head(count)

library(rjags)
library(xtable)
#############modelm - using point mass mixture prior for signals###############
modelm <- "
model{
# likelihood 
for (i in 1:length(y)){
y[i] ~ dpois(lambda[i])
log(lambda[i]) <- alpha[gene[i]] + (-1)^line[i]*tau[gene[i]] + beta[gene[i]]*cov[i] +eps[i]  
}
# prior level 1
for (i in 1:ngene){
alpha[i] ~ dnorm(0,1/sigmaalpha^2)

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
sigmaalpha ~ dunif(0,100)
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
#######Simulation data###############
library(reshape)
sim_data <- function(K, ngene, mualpha, sigmaalpha, 
                     pitau,mutau, sigmatau, pibeta,
                     mubeta, sigmabeta){
  # prior level 1
  x <- rnorm(2*K, 0,1)
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
  dat <- list(y = melt_count$value, 
              gene = melt_count$X1, 
              line = melt_count$line, 
              cov = melt_count$cov, 
              ngene = ngene,
              tau = tau, 
              beta = beta, 
              alpha = alpha)
  return(dat)
  
}


#########function run simulation data #####################3
out <- function(K, ngene, mualpha, sigmaalpha, 
                pitau,mutau, sigmatau, pibeta,
                mubeta, sigmabeta, epstau, epsbeta){
  data <- sim_data(K, ngene, mualpha, sigmaalpha, 
                   pitau,mutau, sigmatau, pibeta,
                   mubeta, sigmabeta)
  
  mm <- jags.model(textConnection(modelm), data[1:5],n.chains = 1)
  resm <- coda.samples(mm, c("tau","alpha","beta",
                             "pitau","pibeta", 
                             "binbeta","bintau", 
                             "sigmatau","sigmabeta"), 2000) 
  
  mm_tau_est <- which(apply(resm[[1]][,paste("bintau[", 1:ngene,"]",sep ="")], 2, 
                            function(x) mean(abs(1-x))) > 0.5)
  
  mm_tau_est_eps <- which(apply(resm[[1]][,paste("tau[", 1:ngene,"]",sep ="")], 2, 
                                function(x) mean(abs(x)>epstau)) > 0.5)
  
  mm_tau_true <- which(data$tau!=0)
  mm_tau_correct <- sum(mm_tau_est%in%mm_tau_true)
  mm_tau_correct_eps <- sum(mm_tau_est_eps%in%mm_tau_true)
  
  mh <- jags.model(textConnection(modelh), data[1:5] ,n.chains = 1)
  resh <- coda.samples(mh, c("tau","alpha","beta",
                             "sigmatauj", "sigmabetaj",  
                             "sigmatau","sigmabeta"), 2000) 
  
  
  mh_tau_est <- which(apply(resh[[1]][,paste("sigmatauj[", 1:ngene,"]",sep ="")], 2, 
                            function(x) mean(1-1/(x^2+1))) > 0.5)
  
  mh_tau_est_eps <- which(apply(resh[[1]][,paste("tau[", 1:ngene,"]",sep ="")], 2, 
                                function(x) mean(abs(x)>epstau)) > 0.5)
  mh_tau_true <- which(data$tau!=0)
  mh_tau_correct <- sum(mh_tau_est%in%mh_tau_true)
  mh_tau_correct_eps <- sum(mh_tau_est_eps%in%mh_tau_true)
  
  return(c(   mm_tau_est = length(mm_tau_est), 
              mm_tau_correct_est = mm_tau_correct,
              mh_tau_est = length(mh_tau_est), 
              mh_tau_correct_est = mh_tau_correct,
              
              mm_tau_est_eps = length(mm_tau_est_eps), 
              mm_tau_correct_est_eps = mm_tau_correct_eps,
              mh_tau_est_eps = length(mh_tau_est_eps), 
              mh_tau_correct_est_eps = mh_tau_correct_eps,
              
              tau_true = length(mm_tau_true)
  ))
}


########### Actually run simulation data ####################
K <- 12

ngene <- 100
# prior level 2
mualpha <- 3
sigmaalpha <- 2
pitau <-  0.8
mutau <- c(0.5,1,2)
sigmatau <- 0.25
pibeta <- 0.8
mubeta <- c(0.5,1,2)
sigmabeta  <- 0.25
epstau <- epsbeta <- 2*mutau/3
post_out <- array(0, dim = c(3,3,9))
for(i in 1:3){
  for (j in 1:3){
    post_out[i,j,] <- out(K, ngene, mualpha, sigmaalpha, 
                          pitau,mutau[i], sigmatau, pibeta,
                          mubeta[j] , sigmabeta, epstau[i], epsbeta[j])
    
  }
}


##########result ########################
library(xtable)
res <- cbind(post_out[1,1,],
             post_out[1,2,],
             post_out[1,3,],
             post_out[2,1,],
             post_out[2,2,],
             post_out[2,3,],
             post_out[3,1,],
             post_out[3,2,],
             post_out[3,3,])
class(res)
colnames(res) <- c("$\\alpha_1$ (0.5,0.5)", "(0.5,1)", "(0.5,2)", 
                   "(1,0.5)", "(1,1)", "(1,2)", 
                   "(2,0.5)", "(2,1)", "(2,2)")


rownames(res) <- c("$\\alpha_1$","mm.c.est.DE.tau",
                   "mh.est.DE.tau","mh.c.est.DE.tau",
                   "mm.est.DE.tau.eps","mm.c.est.DE.tau.eps",
                   "mh.est.DE.tau.eps","mh.c.est.DE.tau.eps",
                   "true.DE.tau")

print(xtable(res, digit = 0, caption = "The meaning of rownames: mm stands for mixture 
             prior, mh stands for horseshoe prior.$\\alpha_1$
             mm.est.DE.tau is the  number of estimated signal tau when using pointmass 
             mixture prior (mm).
             mm.c.est.DE.tau is the number of correctly estimated singal taus when using
             pointmass mixture prior.
             mm.est.DE.tau.eps is the number of estimated signal tau when using pointmass 
             mixture prior and the second metric (varepsilon).
             mm.c.est.DE.tau.eps is the number of correctly estimated 
             signal tau when using pointmass mixture prior and the second metric.
             true.DE.tau is the true number of signal taus (we know this because of
             simulation). 
             "),sanitize.text.function=function(x){x})
@
  
  <<restau,results='hide',echo=FALSE>>=
  library(xtable)
res <- cbind(post_out[1,1,],
             post_out[1,2,],
             post_out[1,3,],
             post_out[2,1,],
             post_out[2,2,],
             post_out[2,3,],
             post_out[3,1,],
             post_out[3,2,],
             post_out[3,3,])
class(res)
colnames(res) <- c("$\\alpha_1$ (0.5,0.5)", "(0.5,1)", "(0.5,2)", 
                   "(1,0.5)", "(1,1)", "(1,2)", 
                   "(2,0.5)", "(2,1)", "(2,2)")


rownames(res) <- c("$\\alpha_1$","mm.c.est.DE.tau",
                   "mh.est.DE.tau","mh.c.est.DE.tau",
                   "mm.est.DE.tau.eps","mm.c.est.DE.tau.eps",
                   "mh.est.DE.tau.eps","mh.c.est.DE.tau.eps",
                   "true.DE.tau")

print(xtable(res, digit = 0, caption = "The meaning of rownames: mm stands for mixture 
       prior, mh stands for horseshoe prior.$\\alpha_1$
       mm.est.DE.tau is the  number of estimated signal tau when using pointmass 
       mixture prior (mm).
       mm.c.est.DE.tau is the number of correctly estimated singal taus when using
       pointmass mixture prior.
       mm.est.DE.tau.eps is the number of estimated signal tau when using pointmass 
       mixture prior and the second metric (varepsilon).
       mm.c.est.DE.tau.eps is the number of correctly estimated 
       signal tau when using pointmass mixture prior and the second metric.
       true.DE.tau is the true number of signal taus (we know this because of
       simulation). 
       "),sanitize.text.function=function(x){x})
