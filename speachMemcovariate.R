source("countfaMCMC-covars.R")

fulldata = read.csv("speachMemcovariateDTM.csv")

# Y: word counts
# X: continuous covariates
Y = fulldata[,15:ncol(fulldata)]

speakernames = fulldata$memberID

conservative = {fulldata[,3] + fulldata[,5]}
X = fulldata[, c(9:12, 14)]
X = cbind(conservative, X)

# Data set properties
N = nrow(Y)
P = ncol(Y)
D = ncol(X)

# How many factors?
K = 10


# Negative-binomial hyperparameter
overdisp = 1

# Load document-term matrix
Y = as.matrix(Y)
X = as.matrix(X)

## Run the MCMC and extract the summaries
ndraws = 200
nburn = 100
thin = 2
mcmc1 = countfactormcmc(Y, X, K=K, nmc=ndraws, burn=nburn, keep=thin, verbose=1)


ahat.counts = colMeans(mcmc1$interceptCounts)
ahat.numerical = colMeans(mcmc1$interceptNumerical)
muhat = colMeans(mcmc1$factormeans) 


## Postprocess each draw
loadings.full = array(0, c(D+P, K, ndraws))
factorscores = array(0, c(N, K, ndraws))

# Set the signs of the factors for identification purposes
thisdraw = rbind(mcmc1$loadingsNumerical[,,1], mcmc1$loadingsCounts[,,1])
thissvd = svd(thisdraw)
ltemp = thissvd$u %*% diag(thissvd$d)
colmax = apply(ltemp,2, max)
colmin = apply(ltemp,2, min)
index.tocheck = rep(0,K)
index.sign = rep(1,K)
for(j in 1:K)
{
	if(abs(colmax[j]) > abs(colmin[j]))
	{
		index.tocheck[j] = which.max(ltemp[,j])
		index.sign[j] = 1
	}
	else
	{
		index.tocheck[j] = which.min(ltemp[,j])
		index.sign[j] = -1
	}
}

for(i in 1:ndraws)
{
	thisdraw = rbind(mcmc1$loadingsNumerical[,,i], mcmc1$loadingsCounts[,,i])
	thissvd = svd(thisdraw)
	ltemp = thissvd$u %*% diag(thissvd$d)
	signcheck = sign(diag(ltemp[index.tocheck,]))
	loadings.full[,,i] = ltemp %*% diag(signcheck * index.sign)
	factorscores[,,i] = mcmc1$factorscores[,,i] %*% {thissvd$v %*% diag(signcheck * index.sign)}
}

loadings.predictors = loadings.postmean[1:D,]
rownames(loadings.predictors) = colnames(X)
loadings.topics = loadings.postmean[(D+1):(D+P),]
rownames(loadings.topics) = colnames(Y)
loadings.predictors
sort(loadings.topics[,1], decreasing=FALSE)[1:20]

factors.postmean = apply(factorscores, c(1,2), mean)



