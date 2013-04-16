library(mvtnorm)
library(BayesLogit)
library(doParallel)
library(foreach)
library(iterators)
library(doMC)
registerDoMC()
library(RcppArmadillo)
sourceCpp("nbtools.cpp")
#registerDoParallel()

NCores = multicore:::detectCores()
	
###### Functions for conditional draws

drawomegamatrix = function(numtrials, linearpredictor)
{

	P = ncol(numtrials)
	chunksize = ceiling(P/NCores)
	# Using a parallel for comprehension from foreach
	omegamatrix = foreach(	nt = iter(numtrials, by='col', chunksize=chunksize),
							eta = iter(linearpredictor, by='col', chunksize=chunksize),
							.combine='cbind') %dopar% {
		drawomegamatrixBlock(nt, eta)
	}
	return(omegamatrix)
}

drawomegamatrixBlock = function(numtrials, linearpredictor)
{
	# Vectorized to C++
	omvec = rpg.devroye(prod(dim(numtrials)), numtrials, as.numeric(linearpredictor))
	return(matrix(omvec, nrow=nrow(numtrials), ncol=ncol(numtrials)))
}


drawloadings.counts = function(Y, factorscores, omegamatrix, factormeans, overdisp, priorprecision=NULL)
# factorscores is an NxK matrix of factor scores
# Y is an NxP matrix of counts, each row a single unit
# omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
{
#	NCores = multicore:::detectCores();
	P = ncol(Y)
	K = ncol(factorscores)
	if(missing(priorprecision)) priorprecision = diag(0.01,K)
	chunksize = ceiling(P/NCores)
	loadings = foreach(	omega = iter(omegamatrix, by='col', chunksize=chunksize),
						y = iter(Y, by='col', chunksize=chunksize),
						.combine='rbind')	%dopar% {		
		#b = matrix(drawloadingsBlockC(y, factorscores, omega, factormeans, overdisp, diag(0,K))$loadings, ncol=K)
		b = drawloadingsBlock.counts(y, factorscores, omega, factormeans, overdisp)
		b
	}
	a = loadings[,1]
	b = loadings[,2:(K+1)]
	list(a=a, b=b)
}

drawloadingsBlock.counts = function(Y, factorscores, omegamatrix, factormeans, overdisp, priorprecision=NULL)
# factorscores is an NxK matrix of factor scores
# Z is an NxP matrix of pseudo-data, each row a single unit
# omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
# This returns the block of the loadings matrix, assuming there is an intercept term (e.g. a factor of one)
{
	P = ncol(Y)
	K = ncol(factorscores)
	design = cbind(1, factorscores)
	if(missing(priorprecision)) priorprecision = diag(1,K+1)
	loadings = matrix(0, nrow=P, ncol=K+1)
	Eye = diag(1,K+1)
	for(j in 1:P)
	{
		omega = omegamatrix[,j]
		kappa = (Y[,j] - overdisp)/2 - omega*factormeans
		#postvar = solve(crossprod(factorscores, omega*factorscores) + priorprecision)
		postprec = crossprod(design, omega*design) + priorprecision
		L = t(chol(postprec))
		Li = forwardsolve(L, Eye)
		postvar = crossprod(Li)
		postmean = postvar %*% {crossprod(design, kappa)}
		loadings[j,] = postmean + crossprod(Li, rnorm(K+1))
	}
	return(loadings)
}

# drawloadingsBlock = function(Y, factorscores, omegamatrix, factormeans, intercept, overdisp, priorprecision=NULL)
# # factorscores is an NxK matrix of factor scores
# # Y is an NxP matrix of observations, each row a single unit
# # omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
# {
	# P = ncol(Y)
	# K = ncol(factorscores)
	# if(missing(priorprecision)) priorprecision = diag(0.01,K)
	# loadings = matrix(0, nrow=P, ncol=K)
	# Eye = diag(1,K)
	# for(j in 1:P)
	# {
		# omega = omegamatrix[,j]
		# kappa = (Y[,j] - overdisp)/2 - omega*(factormeans + intercept[j])
		# #postvar = solve(crossprod(factorscores, omega*factorscores) + priorprecision)
		# postprec = crossprod(factorscores, omega*factorscores) + priorprecision
		# L = t(chol(postprec))
		# Li = forwardsolve(L, Eye)
		# postvar = crossprod(Li)
		# postmean = postvar %*% {crossprod(factorscores, kappa)}
		# loadings[j,] = postmean + crossprod(Li, rnorm(K))
	# }
	# return(loadings)
# }

# drawintercept = function(Y, factorscores, omegamatrix, factormeans, loadings, overdisp)
# # Draw the intercept of length P
# {
	# NCores = multicore:::detectCores()
	# partialpredictor = cbind(loadings,1) %*% rbind(t(factorscores), factormeans)
	# partialpredictor = t(partialpredictor)
	# P = ncol(Y)
	# N = nrow(Y)
	# Z = (Y-overdisp)/2 - partialpredictor*omegamatrix
	# PostPrec = colSums(omegamatrix)
	# PostMeans = colSums(Z)/PostPrec
	# intercept = rnorm(P,PostMeans, 1/sqrt(PostPrec))
	# return(intercept)
# }

drawloadings.numerical = function(X, factorscores, priorprecision=NULL)
{
	D = ncol(X)
	K = ncol(factorscores)
	loadings = matrix(0, nrow=D, ncol=K+1)
	design = cbind(1, factorscores)
	if(missing(priorprecision)) priorprecision = diag(0.01,K+1)
	postprec = crossprod(design) + priorprecision
	postvar = solve(postprec)
	for(j in 1:D)
	{
		postmean = postvar %*% crossprod(design, X[,j])
		loadings[j,] = rmvnorm(1, postmean, postvar)
	}
	a = loadings[,1]
	b = loadings[,2:(K+1)]
	list(a=a,b=b)
}


drawfactormeans = function(Y, factorscores, omegamatrix, intercept, loadings, overdisp)
# Draw the intercept of length P
{
#	NCores = multicore:::detectCores()
	partialpredictor = cbind(intercept,loadings) %*% rbind(1, t(factorscores))
	partialpredictor = t(partialpredictor)
	P = ncol(Y)
	N = nrow(Y)
	
	Z = (Y-overdisp)/2 - partialpredictor*omegamatrix
	PostPrec = rowSums(omegamatrix)
	PostMeans = rowSums(Z)/PostPrec
	factormeans = rnorm(N,PostMeans, 1/sqrt(PostPrec))
	return(factormeans)
}


drawscores = function(Y, X, b.counts, b.numerical, omegamatrix, a.counts, a.numerical, overdisp, factorprecision=NULL)
# b.counts is a PxK matrix of factor loadings for count outcomes
# b.numerical is is a DxK matrix of factor loadings for the numerical outcomes
# Y is an NxP matrix of count observations, each row a single unit
# X is an NxD matrix of numerical observations/covariates, each row a single unit
# omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
# a.counts and a.numerical are the feature-level intercepts
# overdisp is the NB overdispersion parameter
{
#	NCores = multicore:::detectCores()
	N = nrow(Y)
	K = ncol(b.counts)
	if(missing(factorprecision)) factorprecision = diag(1,K+1)
	chunksize = ceiling(N/NCores)
	factorscores = foreach(	om = iter(omegamatrix, by='row', chunksize=chunksize),
						y = iter(Y, by='row', chunksize=chunksize),
						x = iter(X, by='row', chunksize=chunksize),
						.combine='rbind')	%dopar% {
		scores = drawscoresBlock(y, x, b.counts, b.numerical, om, a.counts, a.numerical, overdisp, factorprecision)
		scores;
	}
	factorscores;
	# mu = factorscores[,1]
	# f = factorscores[,2:(K+1)]
	# list(mu=mu, f=f)
}

drawscoresBlock = function(Y, X, b.counts, b.numerical, omegamatrix, a.counts, a.numerical, overdisp, factorprecision=NULL)
# loadings is a PxK matrix of factor loadings
# Y is an NxP matrix of observations, each row a single unit
# omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
{
	N = nrow(Y)
	P = ncol(Y)
	D = ncol(X)
	K = ncol(b.counts)
	Eye = diag(1,K+1)
	if(missing(factorprecision)) factorprecision = diag(1,K+1)
	factorscores = matrix(0, nrow=N, ncol=K+1)
	design = rbind(b.numerical, b.counts)
	design = cbind(1, design)
	for(i in 1:N)
	{
		z = (Y[i,] - overdisp)/{2*omegamatrix[i,]}
		outcome = as.numeric(c(X[i,], z)) - c(a.numerical, a.counts)
		errorprec = as.numeric(c(rep(1, D), omegamatrix[i,]))
		postprec = crossprod(design, errorprec*design) + factorprecision
		L = t(chol(postprec))
		Li = forwardsolve(L, Eye)
		postvar = crossprod(Li)
		postmean = postvar %*% {crossprod(design, errorprec*outcome)}
		factorscores[i,] = postmean + crossprod(Li, rnorm(K+1))
		#factorscores[i,] = rmvnorm(1, postmean, postvar)
	}
	as.matrix(factorscores)
}

# drawscoresBlock = function(Y, loadings, omegamatrix, partialpredictor, overdisp, factorprecision=NULL)
# {
	# N = nrow(Y)
	# K = ncol(loadings)
	# Eye = diag(1,K)
	# if(missing(factorprecision)) factorprecision = diag(1,K)
	# factorscores = matrix(0, nrow=N, ncol=K)
	# for(i in 1:N)
	# {
		# kappa = (Y[i,] - overdisp)/2 - omegamatrix[i,]*partialpredictor[i,]
		# postprec = crossprod(loadings, omegamatrix[i,]*loadings) + factorprecision
		# L = t(chol(postprec))
		# Li = forwardsolve(L, Eye)
		# postvar = crossprod(Li)
		# #postvar = chol2inv(U)
		# #postvar = solve(crossprod(loadings, omega*loadings) + factorprecision)
		# postmean = postvar %*% {crossprod(loadings, kappa)}
		# factorscores[i,] = postmean + crossprod(Li, rnorm(K))
		# #factorscores[i,] = rmvnorm(1, postmean, postvar)
	# }
	# return(factorscores)
# }


### Main MCMC function
countfactormcmc = function(Y, X, K, overdisp=1, burn=50, nmc = 500, keep = 1, verbose=NULL)
# Y is a (numsamples X numfeatures) matrix of count outcomes
# Y is a (numsamples X D) matrix of numerical covariates (right now including binary)
# K is the number of factors
# overdisp is the first parameter in the NB parameterization, ie y ~ NB(overdisp, logit(linearfunc))
{
	if(missing(verbose)) verbose=nmc+burn
	
	# Bookkeeping
	P = ncol(Y)
	D = ncol(X)
	N = nrow(Y)
	numtrials = Y + overdisp
	
	# Set up placeholders to save the relevant draws
	save.loadingsCounts = array(0, dim=c(P,K,nmc))
	save.loadingsNumerical = array(0, dim=c(D,K,nmc))
	save.interceptCounts = matrix(0, nmc, P)
	save.interceptNumerical = matrix(0, nmc, D)
	save.factorscores = array(0, dim=c(N,K,nmc))
	save.factormeans = matrix(0, nmc, N)
	
	# Initialize the chain
	a.counts = rep(0,P)
	b.counts = matrix(rnorm(P*K, 0, 1), nrow=P, ncol=K)
	a.numerical = rep(0,D)
	b.numerical = matrix(rnorm(D*K, 0, 1), nrow=D, ncol=K)
	
	factorscores = matrix(rnorm(N*K,0,1), nrow=N, ncol=K)
	factormeans = rep(0,N)
	
	# Main loop
	for(draw in 1:(nmc*keep+burn))
	{
		if(draw %% verbose == 0) cat("\tDraw", draw, "\n")
		
		linearpredictor = cbind(a.counts,b.counts,1) %*% rbind(1,t(factorscores), factormeans)
		linearpredictor = t(linearpredictor)
		
		omega = drawomegamatrix(numtrials, linearpredictor)
		loadings.counts = drawloadings.counts(Y, factorscores, omega, factormeans, overdisp)
		a.counts = loadings.counts$a; b.counts = loadings.counts$b
		
		loadings.numerical = drawloadings.numerical(X, factorscores)
		a.numerical = loadings.numerical$a; b.numerical = loadings.numerical$b
		
		scoredraw = drawscores(Y, X, b.counts, b.numerical, omega, a.counts, a.numerical, overdisp)
		factormeans = scoredraw[,1]
		factorscores = scoredraw[,2:(K+1)]

		# save draws
		if(draw > burn && (draw - burn) %% keep == 0)
		{
			save.loadingsCounts[,,(draw-burn)/keep] = b.counts
			save.loadingsNumerical[,,(draw-burn)/keep] = b.numerical
			save.interceptCounts[(draw-burn)/keep,] = a.counts
			save.interceptNumerical[(draw-burn)/keep,] = a.numerical
			save.factorscores[,,(draw-burn)/keep] = factorscores
			save.factormeans[(draw-burn)/keep,] = factormeans
		}
	}
	list(	loadingsCounts = save.loadingsCounts, loadingsNumerical = save.loadingsNumerical,
			interceptCounts = save.interceptCounts, interceptNumerical = save.interceptNumerical,
			factorscores = save.factorscores, factormeans = save.factormeans)
}
