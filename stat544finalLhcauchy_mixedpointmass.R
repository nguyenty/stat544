
library(rjags)

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

  tau[j] ~ dnorm(0, 1/st[j]^2)
  st[j] ~ dt(0, 1/sct^2, 1) T(0,)
  
  beta[j] ~ dnorm(0, 1/sb[j]^2)
  sb[j] ~ dt(0, 1/scb^2, 1) T(0,)
}
# prior level 2
sct ~ dt(0, 1, 1) T(0,)
scb ~ dt(0, 1, 1) T(0,)
}
"
#load(file = "/home/ntyet/New folder/STAT544/HW/finalproject/data.RData")

#load(file = "/run/user/1000/gvfs/smb-share:server=smb.stat.iastate.edu,share=ntyet/New folder/STAT544/HW/finalproject/data.RData")


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
# prior level 2
pitau ~ dbeta(8,1)
pibeta ~ dbeta(8,1)
sigmatau ~ dunif(0,100)
sigmabeta ~ dunif(0,100)
}
"

# getwd()
# save(data, file = "data.RData")
# load("data.RData")
# data
# unique(gene)
# unique(line)
# unique(cov)
# length(cov)
# length(gene)
# length(line)
load(file = "P:/research/K_20/pbeta_0.5/p.beta_0.5i.beta_0.1e.beta_0.5m_1.RData")
str(sim1)
counts <- sim1$counts
x <- sim1$x

#first 200 genes
ngene <- dim(counts)[1]
m <- dim(counts)[2]
count <- counts[1:ngene,]
head(count)
colnames(count)
rownames(count) <- as.character(1:ngene)
head(melt(count))
Line <- rep(c(1,2), each = m/2)
colnames(count)
colnames(count)[Line ==1]
RFI
unique(melt(count)$X2)
melt_count <- melt(count)
melt_count$line <- 1
head(melt_count)
melt_count$line[melt_count$X2%in%Line[(m/2+1):m]] <- 2
melt_count$cov <- 1
for(i in 1:m) melt_count$cov[melt_count$X2==i] <- x[i]
# head(melt_count)
# RFI
dat <- as.data.frame(cbind(gene = melt_count$X1, 
                           line = melt_count$line, 
                           cov = melt_count$cov, 
                           value = melt_count$value))

head(dat)
data <- list(y = dat$value, 
             gene = dat$gene, 
             line = dat$line, 
             cov = dat$cov, 
             ngene = ngene)

m0 <- proc.time()
mm <- jags.model(textConnection(model2), data, n.adapt=1000,n.chains = 1) # mix point mass

resm <- coda.samples(mm, c("tau", "beta"), 5000) # mix point mass
proc.time() -m0

m0 <- proc.time()
mh <- jags.model(textConnection(model1), data, n.chains = 1) # horseshoe

resh <- coda.samples(mh, c("tau", "beta"), 5000) # horseshoe
proc.time() -m0

# #summary(res)
# str(res[[1]])
# str(res)
# count[2,Line ==1]
# count[2,Line ==2]
# summary(res[[1]][,2])
# plot(res)

plot(res[, c("tau[1]", "tau[2]", 
             "tau[5]", "tau[1000]")])
mean(abs(res[[1]][,"tau[993]"])> .1)
#sim1$beta.ind
plot(res1[, c("tau[1]", "tau[2]", 
             "tau[3]", "tau[4]")])
mean(abs(res1[[1]][,"tau[100]"])< .1)
(names(apply(res1[[1]], 2, mean)))
sum(apply(res1[[1]], 2, mean)[grep("tau", names(apply(res1[[1]], 2, mean)))]<0.1)

sum(apply(res[[1]], 2, mean)[grep("tau", names(apply(res[[1]], 2, mean)))]<0.1)
which((apply(res[[1]], 2, mean)[grep("tau", names(apply(res[[1]], 2, mean)))]>0.2))

count[20,Line ==1]
count[20,Line ==2]


##########check using QuasiSeq
K <- 20
design.list <- vector("list", 3)
x1 <- as.factor(rep(1:2, each = K))

design.list[[1]] <- model.matrix(~ x1 + x)
design.list[[2]] <- model.matrix(~ x1) # test for covariate
design.list[[3]] <- model.matrix(~ x) # test for treatment effect
test.mat <- rbind(1:2, c(1,3))
row.names(test.mat) <- c("Covariate", "Treatment")
fit.2 <- QL.fit(counts, design.list, 
                test.mat,
                Model="NegBin", 
                print.progress=FALSE)

result.fit2 <- QL.results(fit.2,Plot= FALSE)

pvalue.cov <- result.fit2$P.values[[3]][,1]
pvalue.trt.cov <- result.fit2$P.values[[3]][,2]
hist(pvalue.cov, nclass = 200)
hist(pvalue.trt.cov, nclass = 200)
str(result.fit2)
result.fit2$m0
apply(result.fit2$Q.values[[3]], 2, function(x) which(x<=0.05))

post_tau <- apply(res[[1]][,1001:2000],2, function(x) mean(x>0.1)) # mix point mass
sum((post_tau >0.5))
sum(which(post_tau >0.5)>200)
count[1000,]

post_tau1 <- apply(res1[[1]][,1001:2000],2, function(x) mean(x>0.01)) # horseshoe prior
sum((post_tau1 >0.1))
sum(which(post_tau1 >0.5)>200)


# simulated data from  point mass mixture model 
# covariate 
library(reshape)
K <- 10
x <- rnorm(2*K, 0,1)
ngene <- 100
# prior level 2
pitau <-  0.8
pibeta <- 0.8
sigmatau <- 1
sigmabeta  <- 1

# prior level 1

tau <-(1-rbinom(ngene,1,pitau))*rnorm(ngene, 4, sigmatau)
beta <- (1-rbinom(ngene,1,pibeta))*rnorm(ngene, 1, sigmabeta)
alpha <- rnorm(ngene, 2,1 )
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

head(count)
head(lambda)
count[1,]
lambda[1,]

rownames(count) <- as.character(1:ngene)
head(melt(count))
line <- rep(c(1,2), each = K)

unique(melt(count)$X2)
melt_count <- melt(count)
melt_count$line <- 1
head(melt_count)
melt_count$line[melt_count$X2%in%line[(K+1):(2*K)]] <- 2
melt_count$cov <- 1
for(i in 1:(2*K)) melt_count$cov[melt_count$X2==i] <- x[i]

head(melt_count)
count[1:10,]
data <- list(y = melt_count$value, 
             gene = melt_count$X1, 
             line = melt_count$line, 
             cov = melt_count$cov, 
             ngene = ngene)

m0 <- proc.time()
mm <- jags.model(textConnection(modelm), data,n.chains = 1) # mix point mass
#?jags.model
resm <- coda.samples(mm, c("tau",  "alpha","beta"), 5000) # mix point mass
proc.time() -m0

plot(resm[, c("tau[25]", "tau[21]")])
which(tau!=0)
tau[which(tau!=0)]
alpha
plot(resm[, c("alpha[1]", "alpha[2]")])

plot(resm[, c("beta[4]", "beta[49]")])
which(beta!=0)
beta[which(beta!=0)]
beta[4]
m0 <- proc.time()
mh <- jags.model(textConnection(modelh), data, n.chains = 1) # horseshoe

resh <- coda.samples(mh, c("tau", "beta", "alpha"), 5000) # horseshoe
proc.time() -m0

plot(resh[, c("tau[3]", "tau[2]")])
which(tau!=0)
tau[which(tau!=0)]
alpha
plot(resm[, c("alpha[1]", "alpha[2]")])

plot(resm[, c("beta[2]", "beta[49]")])
which(beta!=0)
beta[which(beta!=0)]


#sim1$beta.ind
plot(resh[, c("tau[46]", "tau[47]")])
tau
which(tau!=0)
post_taum <- apply(resm[[1]][,51:100],2, function(x) mean(x>0.05)) # mix point mass
sum((post_taum >0.5))
sum(which(post_tau >0.5)>200)
count[1000,]
post_taum[which(tau!=0)]

post_tauh <- apply(resh[[1]][,51:100],2, function(x) mean(x>0.01)) # horseshoe prior
sum((post_tauh >0.1))
sum(which(post_tau1 >0.5)>200)