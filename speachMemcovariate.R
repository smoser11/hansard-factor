source("countfaMCMC-covars.R")

fulldata = read.csv("speachMemcovariateDTM.csv")

# Y: word counts
# X: continuous covariates
Y = fulldata[,15:ncol(fulldata)]
X = fulldata[, c(3:6, 9:12, 14)]

# Data set properties
N = fulldata$nrow
P = fulldata$ncol

# How many factors?
K = 10


# Negative-binomial hyperparameter
overdisp = 1

# Load document-term matrix
Y = as.matrix(Y)
X = as.matrix(X)

## Run the MCMC and extract the summaries
mcmc1 = countfactormcmc(Y, X, K=K, nmc=10, burn=0, verbose=1)

bhat = apply(mcmc1$loadings,c(1,2), mean)
ahat = colMeans(mcmc1$intercept)
muhat = colMeans(mcmc1$factormeans) 

