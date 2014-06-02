library(rjags)
library(xtable)
#############modelm - using point mass mixture prior for signals###############
modelm <- "
model{
# likelihood 
for (i in 1:length(y)){
y[i] ~ dnegbin((1/omega[gene[i]])/(lambda[i] + 1/omega[gene[i]]), 1/omega[gene[i]])
log(lambda[i]) <- alpha[gene[i]] + (-1)^line[i]*tau[gene[i]] + beta[gene[i]]*cov[i]   
}
# prior level 1
for (i in 1:ngene){
alpha[i] ~ dnorm(0,1/10000)
omega[i] ~ dlnorm(0, 1/10000)
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



#############Sim_data_function ###################
library(reshape)
sim_data <- function(K, ngene, mualpha, sigmaalpha, 
                     pitau,mutau, sigmatau, pibeta,
                     mubeta, sigmabeta){
  # prior level 1 
  x <- rnorm(2*K, 0,1)
  bintau <- rbinom(ngene,1,pitau)
  tau <-(1-bintau)*rnorm(ngene, mutau, sigmatau)
  binbeta <- rbinom(ngene,1,pibeta)
  beta <- (1-binbeta)*rnorm(ngene, mubeta, sigmabeta)
  alpha <- rnorm(ngene, mualpha,sigmaalpha )
  lambda <- matrix(0, ncol = 2*K, nrow = ngene)
  omega <- exp(rnorm(ngene, 0, 2))
  count <- matrix(0, ncol = 2*K, nrow = ngene)
  
  for (j in 1:ngene){
    for (k in 1:K){
      lambda[j,k] <- exp(alpha[j] - tau[j] + beta[j]*x[k])
      lambda[j, k+K] <- exp(alpha[j] + tau[j] + beta[j]*x[k+K])
      count[j,k] <- rnbinom(1, size = 1/omega[j], mu = lambda[j,k])
      count[j,k+K] <- rnbinom(1, size = 1/omega[j], mu = lambda[j,k+K])
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
              bintau=bintau,
              tau = tau, 
              binbeta = binbeta,
              beta = beta, 
              alpha = alpha, 
              omega = omega)
  return(dat)
  
}


###############run_mm_simulationdata#######################
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
    
#    mh <- jags.model(textConnection(modelh), data[1:5] ,n.chains = 1)
#    resh <- coda.samples(mh, c("tau","alpha","beta",
#                               "sigmatauj", "sigmabetaj",  
#                               "sigmatau","sigmabeta"), 2000) 
    
    
#    mh_tau_est <- which(apply(resh[[1]][,paste("sigmatauj[", 1:ngene,"]",sep ="")], 2, 
#                              function(x) mean(1-1/(x^2+1))) > 0.5)
    
#    mh_tau_est_eps <- which(apply(resh[[1]][,paste("tau[", 1:ngene,"]",sep ="")], 2, 
#                                  function(x) mean(abs(x)>epstau)) > 0.5)
#    mh_tau_true <- which(data$tau!=0)
#    mh_tau_correct <- sum(mh_tau_est%in%mh_tau_true)
#    mh_tau_correct_eps <- sum(mh_tau_est_eps%in%mh_tau_true)
    
    return(c(   mm_tau_est = length(mm_tau_est), 
                mm_tau_correct_est = mm_tau_correct,
#                mh_tau_est = length(mh_tau_est), 
 #               mh_tau_correct_est = mh_tau_correct,
                
                mm_tau_est_eps = length(mm_tau_est_eps), 
                mm_tau_correct_est_eps = mm_tau_correct_eps,
  #              mh_tau_est_eps = length(mh_tau_est_eps), 
   #             mh_tau_correct_est_eps = mh_tau_correct_eps,
                
                tau_true = length(mm_tau_true)
    ))
  }

@
  
  
######run_sim##########
K <- 12
ngene <- 100
# prior level 2
mualpha <- 3
sigmaalpha <- 2
pitau <-  0.8
#mutau <- c(0.5,1,2)
mutau <- 1
sigmatau <- 0.25
pibeta <- 0.8
#mubeta <- c(0.5,1,2)
mubeta <- 1
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
