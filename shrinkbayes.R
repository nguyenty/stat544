library(ShrinkBayes)
load("data.RData")
dat2 <- as.data.frame(cbind(data[[1]], data[[2]], data[[3]], data[[4]]))


colnames(dat2) <- c("y", "gene", "line", "cov")


library(reshape)
dat3 <- cast(dat2, formula = gene ~ cov,  value= "y")
str(dat3)
dim(dat3)
colnames(dat3)
head(dat2)
tail(dat2)
count <- dat3[,-1]
dim(count)
head(count)
cov <- as.numeric(colnames(count))
line <- rep(c(1,2), each=12)
line
head(count)
count[26,]
count[81,]
counts <- count
head(count)
