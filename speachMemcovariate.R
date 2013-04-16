source("countfaMCMC-covars.R")

fulldata = read.csv("speachMemcovariateDTM.csv")  #2310 speaches by 8468 words + covariates (in columns c(2:12, 8480) )

# Y: word counts, 
# X: continuous covariates

Y = fulldata[,13:ncol(fulldata)-1] #last col of fulldata is 0/1 year67 (see below)

speakernames = fulldata$memberID
#4=party2c (1 for Conservative, 0 for Liberal)
#5:6 = "type.xcounty" and  "type.xuniversity" for constituency type
#7:9 = "country.xireland"  "country.xscotland" "country.xwales"  for region
#10=malaportionment (>0 means overrepresented in HoC, <0 means underrepresented)
#11=miniTRUE (1 if speaker was a minsiter, 0 else)
#12= popchg (percent increase or decrease of population from 1851 to 1861)
#8480 = year67 (1 if speech was given in 1867, 0 if given in 1866)
X = fulldata[, c(4,5:6,10,12)] 
X$year67<-fulldata[,8480]

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
nburn = 50
thin = 5
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
loadings.postmean = apply(loadings.full, c(1,2), mean)
factors.postmean = apply(factorscores,c(1,2),mean)


loadings.predictors = loadings.postmean[1:D,]
rownames(loadings.predictors) = colnames(X)
loadings.topics = loadings.postmean[(D+1):(D+P),]
rownames(loadings.topics) = colnames(Y)

loadings.predictors*100

# Extract the topics
L = 20
TopicsTopL = data.frame(Word = 1:L)
TopicsBottomL = data.frame(Word = 1:L)
for(k in 1:K)
{
	myorder = order(loadings.topics[,k], decreasing=TRUE)
	thiscommand1 = paste0("TopicsTopL$Topic",k," = colnames(Y)[head(myorder,n=L)]")
	eval(parse(text=thiscommand1))
	thiscommand2 = paste0("TopicsBottomL$Topic",k," = colnames(Y)[rev(tail(myorder,n=L))]")
	eval(parse(text=thiscommand2))
}


# Compute locations of speeches in factor space
# Compare the ratio of determinants of the posterior covariance matrix of the factors
# This is a measure of the volume of factor scores K-space for labour
# versus those of conservatives
LtoCVolumeRatio = rep(0, ndraws)
consind = which(X[,1]==1)
for(i in 1:ndraws)
{
	speeches.cons = factorscores[consind,,i]
	speeches.lab = factorscores[-consind,,i]
	LtoCVolumeRatio[i] = det(cov(speeches.lab))/det(cov(speeches.cons))
}
hist(LtoCVolumeRatio)


# Plot locations of speeches in two of the K dimensions
i = 1; j=9
plot(factors.postmean[,c(i,j)])
points(factors.postmean[consind,c(i,j)], col='red', pch=19)
points(factors.postmean[-consind,c(i,j)], col='blue', pch=19)

# Pick out the most and least Conservative words by looking
# at the posterior mean covariance matrix of the latent linear predictors
PostCovHat = tcrossprod(loadings.postmean)
conswordorder = order(PostCovHat[1,(D+1):(D+P)], decreasing=TRUE)
head(colnames(Y)[conswordorder], 20)
tail(colnames(Y)[conswordorder], 20)
