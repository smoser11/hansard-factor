library(mvtnorm)
library(BayesLogit)
library(doParallel)
library(foreach)
library(iterators)
#library(doMC)
#registerDoMC()
registerDoParallel()

NCores = multicore:::detectCores()
NCores = 2
library(mvtnorm)
library(BayesLogit)
library(doParallel)
library(foreach)
library(iterators)
#library(doMC)
#registerDoMC()
registerDoParallel()

NCores = multicore:::detectCores()
NCores = 2

###### Functions for conditional draws

drawthiscol = function(j, nt, lp)
{
  rpg.devroye(nrow(nt), nt[,j], lp[,j])
}

drawomegamatrix = function(numtrials, linearpredictor)
{
  #NCores = multicore:::detectCores()
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

drawloadings = function(Y, factorscores, omegamatrix, factormeans, intercept, overdisp, priorprecision=NULL)
  # factorscores is an NxK matrix of factor scores
  # Y is an NxP matrix of observations, each row a single unit
  # omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
{
  #	NCores = multicore:::detectCores()
  P = ncol(Y)
  K = ncol(factorscores)
  if(missing(priorprecision)) priorprecision = diag(0.01,K)
  
  chunksize = ceiling(P/NCores)
  alpha = as.matrix(intercept)
  loadings = foreach(	omega = iter(omegamatrix, by='col', chunksize=chunksize),
                      y = iter(Y, by='col', chunksize=chunksize),
                      a = iter(alpha,by='row',chunksize=chunksize),
                      .combine='rbind')	%dopar% {		
                        #b = matrix(drawloadingsBlockC(y, factorscores, omega, factormeans, a, overdisp, diag(0,K))$loadings, ncol=K)
                        b = drawloadingsBlock(y, factorscores, omega, factormeans, a, overdisp)
                        b
                      }
  return(loadings)
}


drawloadingsBlock = function(Y, factorscores, omegamatrix, factormeans, intercept, overdisp, priorprecision=NULL)
  # factorscores is an NxK matrix of factor scores
  # Y is an NxP matrix of observations, each row a single unit
  # omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
{
  P = ncol(Y)
  K = ncol(factorscores)
  if(missing(priorprecision)) priorprecision = diag(0.01,K)
  loadings = matrix(0, nrow=P, ncol=K)
  Eye = diag(1,K)
  for(j in 1:P)
  {
    omega = omegamatrix[,j]
    kappa = (Y[,j] - overdisp)/2 - omega*(factormeans + intercept[j])
    #postvar = solve(crossprod(factorscores, omega*factorscores) + priorprecision)
    postprec = crossprod(factorscores, omega*factorscores) + priorprecision
    L = t(chol(postprec))
    Li = forwardsolve(L, Eye)
    postvar = crossprod(Li)
    postmean = postvar %*% {crossprod(factorscores, kappa)}
    loadings[j,] = postmean + crossprod(Li, rnorm(K))
  }
  return(loadings)
}




drawintercept = function(Y, factorscores, omegamatrix, factormeans, loadings, overdisp)
  # Draw the intercept of length P
{
  #	NCores = multicore:::detectCores()
  partialpredictor = cbind(loadings,1) %*% rbind(t(factorscores), factormeans)
  partialpredictor = t(partialpredictor)
  P = ncol(Y)
  N = nrow(Y)
  
  Z = (Y-overdisp)/2 - partialpredictor*omegamatrix
  PostPrec = colSums(omegamatrix)
  PostMeans = colSums(Z)/PostPrec
  intercept = rnorm(P,PostMeans, 1/sqrt(PostPrec))
  
  return(intercept)
  
  # Splitting and sending to separate cores is slower!
  # chunksize = ceiling(P/NCores)
  # intercept = foreach(	y = iter(Y, by='col', chunksize=chunksize),
  # omega = iter(omegamatrix, by='col', chunksize=chunksize),
  # eta = iter(partialpredictor, by='col', chunksize=chunksize),
  # .combine='c') %dopar% {
  # drawinterceptBlock(y, omega, eta, overdisp)
  # }
}

# drawinterceptBlock = function(Y, omegamatrix, partialpredictor, overdisp)
# # Draw the intercept of length P
# {
# P = ncol(Y)
# N = nrow(Y)
# intercept = rep(0,P)
# for(j in 1:P)
# {
# omega = omegamatrix[,j]
# z = (Y[,j] - overdisp)/{2*omega} - partialpredictor[,j]
# postvar = 1/sum(omega)
# postmean = weighted.mean(z, w=omega)
# intercept[j] = rnorm(1, postmean, sqrt(postvar))
# }
# return(intercept)
# }


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
  
  # The naive looping way
  # factormeans = rep(0,N)
  # for(i in 1:N)
  # {
  # omega = omegamatrix[i,]
  # z = (Y[i,] - overdisp)/{2*omega} - partialpredictor[i,]
  # postvar = 1/sum(omega)
  # postmean = weighted.mean(z, w=omega)
  # factormeans[i] = rnorm(1, postmean, sqrt(postvar))
  # }
  
  # Splitting and sending to separate cores is slower!
  # chunksize = ceiling(N/NCores)
  # factormeans = foreach(	y = iter(Y, by='row', chunksize=chunksize),
  # omega = iter(omegamatrix, by='row', chunksize=chunksize),
  # eta = iter(partialpredictor, by='row', chunksize=chunksize),
  # .combine='c') %dopar% {
  # drawfactormeansBlock(y, omega, eta, overdisp)
  # }
  
  return(factormeans)
}


# drawfactormeansBlock = function(Y, omegamatrix, partialpredictor, overdisp)
# # Draw the observation-specific means of length N
# {
# N = nrow(Y)
# factormeans = rep(0,N)
# for(i in 1:N)
# {
# omega = omegamatrix[i,]
# z = (Y[i,] - overdisp)/{2*omega} - partialpredictor[i,]
# postvar = 1/sum(omega)
# postmean = weighted.mean(z, w=omega)
# factormeans[i] = rnorm(1, postmean, sqrt(postvar))
# }

# # Z = (Y-overdisp)/2 - partialpredictor*omegamatrix
# # PostPrec = rowSums(omegamatrix)
# # PostMeans = rowSums(Z)/PostPrec
# # factormeans = rnorm(N,PostMeans, 1/sqrt(PostPrec))
# return(factormeans)
# }


drawscores = function(Y, loadings, omegamatrix, factormeans, intercept, overdisp, qs, nuis, factorprecision=NULL)
  # loadings is a PxK matrix of factor loadings
  # Y is an NxP matrix of observations, each row a single unit
  # omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
{
  #	NCores = multicore:::detectCores()
  N = nrow(Y)
  K = ncol(loadings)
  if(missing(factorprecision)) factorprecision = diag(1,K)
  partialpredictor = cbind(intercept,1) %*% rbind(1, factormeans)
  partialpredictor = t(partialpredictor)
  chunksize = ceiling(N/NCores)
  factorscores = foreach(omega = iter(omegamatrix, by='row', chunksize=chunksize),
                         y = iter(Y, by='row', chunksize=chunksize),
                         eta = iter(partialpredictor,by='row',chunksize=chunksize),
                         .combine='rbind')	%dopar% {
                           drawscoresBlock(y, loadings, omega, eta, overdisp,  qs, nuis,factorprecision)
                         }
  return(factorscores)
}


drawscoresBlock = function(Y, loadings, omegamatrix, partialpredictor, overdisp, qs, nuis, factorprecision=NULL)
{
  N = nrow(Y)
  K = ncol(loadings)
  Eye = diag(1,K)
  if(missing(factorprecision)) factorprecision = diag(1,K)
  factorscores = matrix(0, nrow=N, ncol=K)
  for(i in 1:N)
  {
    kappa = (Y[i,] - overdisp)/2 - omegamatrix[i,]*partialpredictor[i,]
    postprec = crossprod(loadings, omegamatrix[i,]*loadings) + factorprecision
    L = t(chol(postprec))
    Li = forwardsolve(L, Eye)
    postvar = crossprod(Li)
    #postvar = chol2inv(U)
    #postvar = solve(crossprod(loadings, omega*loadings) + factorprecision)
    postmean = postvar %*% {crossprod(loadings, kappa)}
    #sparcity here
    prob1 = diag(dnorm(0,postmean,postvar))
    p1<-prob1 > 0.00001
    qhat<-rep(1,K)
    rhats<-dnorm(0,0,nuis)/prob1
    qhat[p1] <- rhats[p1]/((1-qs[p1])/qs[p1]+rhats[p1])
    
    factorscores[i,] = ifelse(runif(K)<qhat,postmean + crossprod(Li, rnorm(K)),0)
    #factorscores[i,] = rmvnorm(1, postmean, postvar)
  }
  return(factorscores)
}



### Main MCMC function
countfactormcmc = function(Y, K, overdisp=1, burn=50, nmc = 500, qs=NULL, verbose=NULL)
  # Y is a (numsamples X numfeatures) matrix of count outcomes
  # K is the number of factors
  # overdisp is the first parameter in the NB parameterization, ie y ~ NB(overdisp, logit(linearfunc))
{
  if(missing(verbose)) verbose=nmc+burn
  
  # Bookkeeping
  P = ncol(Y)
  N = nrow(Y)
  numtrials = Y + overdisp
  
  nuis<-rep(0,K)
  qs<-rep(0,K)
  # Set up placeholders to save the relevant draws
  save.loadings = array(0, dim=c(P,K,nmc))
  save.intercept = matrix(0, nmc, P)
  save.factormeans = matrix(0, nmc, N)
  save.factorscores = array(0, dim=c(N,K,nmc))
  
  # Initialize the chain
  loadings = matrix(rnorm(P*K, 0, 1), nrow=P, ncol=K)
  factorscores = matrix(rnorm(N*K,0,1), nrow=N, ncol=K)
  intercept = rep(0,P)
  factormeans = rep(0,N)
  
  # Main loop
  for(draw in 1:(nmc+burn))
  {
    if(draw %% verbose == 0) cat("\tDraw", draw, "\n")
    
    linearpredictor = cbind(intercept,loadings,1) %*% rbind(1,t(factorscores), factormeans)
    linearpredictor = t(linearpredictor)
    
    omega = drawomegamatrix(numtrials, linearpredictor)
    loadings = drawloadings(Y, factorscores, omega, factormeans, intercept, overdisp)
    intercept = drawintercept(Y, factorscores, omega, factormeans, loadings, overdisp)
    factormeans = drawfactormeans(Y, factorscores, omega, intercept, loadings, overdisp)
    
    ms<- apply(factorscores,2,function(x) {sum(abs(x)>.0000001)})
    #sjm get q_s for each column
    for(s in 1:K){
      nuis[s]<-  1/(rgamma(1,(2+ms[s])/2,(2+as.numeric(t(factorscores[,s])%*%factorscores[,s]))/2))      
      qs[s]<- rbeta(1,c_s+ms[s],c_s*d_s/2+P-s-ms[s])
    }
    
    factorscores = drawscores(Y, loadings, omega, factormeans, intercept, overdisp, qs, nuis)
    
    # save draws
    if(draw > burn)
    {
      save.factorscores[,,draw-burn] = factorscores
      save.loadings[,,draw-burn] = loadings
      save.intercept[draw-burn,] = intercept
      save.factormeans[draw-burn,] = factormeans
    }
  }
  list(loadings = save.loadings, intercept = save.intercept, factormeans = save.factormeans, factorscores=save.factorscores)
}

###### Functions for conditional draws

drawthiscol = function(j, nt, lp)
{
	rpg.devroye(nrow(nt), nt[,j], lp[,j])
}

drawomegamatrix = function(numtrials, linearpredictor)
{
	#NCores = multicore:::detectCores()
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

drawloadings = function(Y, factorscores, omegamatrix, factormeans, intercept, overdisp, priorprecision=NULL)
# factorscores is an NxK matrix of factor scores
# Y is an NxP matrix of observations, each row a single unit
# omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
{
#	NCores = multicore:::detectCores()
	P = ncol(Y)
	K = ncol(factorscores)
	if(missing(priorprecision)) priorprecision = diag(0.01,K)
	
	chunksize = ceiling(P/NCores)
	alpha = as.matrix(intercept)
	loadings = foreach(	omega = iter(omegamatrix, by='col', chunksize=chunksize),
						y = iter(Y, by='col', chunksize=chunksize),
						a = iter(alpha,by='row',chunksize=chunksize),
						.combine='rbind')	%dopar% {		
		#b = matrix(drawloadingsBlockC(y, factorscores, omega, factormeans, a, overdisp, diag(0,K))$loadings, ncol=K)
		b = drawloadingsBlock(y, factorscores, omega, factormeans, a, overdisp)
		b
	}
	return(loadings)
}


drawloadingsBlock = function(Y, factorscores, omegamatrix, factormeans, intercept, overdisp, priorprecision=NULL)
# factorscores is an NxK matrix of factor scores
# Y is an NxP matrix of observations, each row a single unit
# omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
{
	P = ncol(Y)
	K = ncol(factorscores)
	if(missing(priorprecision)) priorprecision = diag(0.01,K)
	loadings = matrix(0, nrow=P, ncol=K)
	Eye = diag(1,K)
	for(j in 1:P)
	{
		omega = omegamatrix[,j]
		kappa = (Y[,j] - overdisp)/2 - omega*(factormeans + intercept[j])
		#postvar = solve(crossprod(factorscores, omega*factorscores) + priorprecision)
		postprec = crossprod(factorscores, omega*factorscores) + priorprecision
		L = t(chol(postprec))
		Li = forwardsolve(L, Eye)
		postvar = crossprod(Li)
		postmean = postvar %*% {crossprod(factorscores, kappa)}
		loadings[j,] = postmean + crossprod(Li, rnorm(K))
	}
	return(loadings)
}




drawintercept = function(Y, factorscores, omegamatrix, factormeans, loadings, overdisp)
# Draw the intercept of length P
{
#	NCores = multicore:::detectCores()
	partialpredictor = cbind(loadings,1) %*% rbind(t(factorscores), factormeans)
	partialpredictor = t(partialpredictor)
	P = ncol(Y)
	N = nrow(Y)
	
	Z = (Y-overdisp)/2 - partialpredictor*omegamatrix
	PostPrec = colSums(omegamatrix)
	PostMeans = colSums(Z)/PostPrec
	intercept = rnorm(P,PostMeans, 1/sqrt(PostPrec))

	return(intercept)
	
	# Splitting and sending to separate cores is slower!
	# chunksize = ceiling(P/NCores)
	# intercept = foreach(	y = iter(Y, by='col', chunksize=chunksize),
							# omega = iter(omegamatrix, by='col', chunksize=chunksize),
							# eta = iter(partialpredictor, by='col', chunksize=chunksize),
							# .combine='c') %dopar% {
		# drawinterceptBlock(y, omega, eta, overdisp)
	# }
}

# drawinterceptBlock = function(Y, omegamatrix, partialpredictor, overdisp)
# # Draw the intercept of length P
# {
	# P = ncol(Y)
	# N = nrow(Y)
	# intercept = rep(0,P)
	# for(j in 1:P)
	# {
		# omega = omegamatrix[,j]
		# z = (Y[,j] - overdisp)/{2*omega} - partialpredictor[,j]
		# postvar = 1/sum(omega)
		# postmean = weighted.mean(z, w=omega)
		# intercept[j] = rnorm(1, postmean, sqrt(postvar))
	# }
	# return(intercept)
# }


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
	
	# The naive looping way
	# factormeans = rep(0,N)
	# for(i in 1:N)
	# {
		# omega = omegamatrix[i,]
		# z = (Y[i,] - overdisp)/{2*omega} - partialpredictor[i,]
		# postvar = 1/sum(omega)
		# postmean = weighted.mean(z, w=omega)
		# factormeans[i] = rnorm(1, postmean, sqrt(postvar))
	# }
	
	# Splitting and sending to separate cores is slower!
	# chunksize = ceiling(N/NCores)
	# factormeans = foreach(	y = iter(Y, by='row', chunksize=chunksize),
							# omega = iter(omegamatrix, by='row', chunksize=chunksize),
							# eta = iter(partialpredictor, by='row', chunksize=chunksize),
							# .combine='c') %dopar% {
		# drawfactormeansBlock(y, omega, eta, overdisp)
	# }
	
	return(factormeans)
}


# drawfactormeansBlock = function(Y, omegamatrix, partialpredictor, overdisp)
# # Draw the observation-specific means of length N
# {
	# N = nrow(Y)
	# factormeans = rep(0,N)
	# for(i in 1:N)
	# {
		# omega = omegamatrix[i,]
		# z = (Y[i,] - overdisp)/{2*omega} - partialpredictor[i,]
		# postvar = 1/sum(omega)
		# postmean = weighted.mean(z, w=omega)
		# factormeans[i] = rnorm(1, postmean, sqrt(postvar))
	# }
	
	# # Z = (Y-overdisp)/2 - partialpredictor*omegamatrix
	# # PostPrec = rowSums(omegamatrix)
	# # PostMeans = rowSums(Z)/PostPrec
	# # factormeans = rnorm(N,PostMeans, 1/sqrt(PostPrec))
	# return(factormeans)
# }


drawscores = function(Y, loadings, omegamatrix, factormeans, intercept, overdisp, qs, nuis, factorprecision=NULL)
# loadings is a PxK matrix of factor loadings
# Y is an NxP matrix of observations, each row a single unit
# omegamatrix is an NxP matrix of augmentation variables, one-to-one with Y
{
#	NCores = multicore:::detectCores()
	N = nrow(Y)
	K = ncol(loadings)
	if(missing(factorprecision)) factorprecision = diag(1,K)
	partialpredictor = cbind(intercept,1) %*% rbind(1, factormeans)
	partialpredictor = t(partialpredictor)
	chunksize = ceiling(N/NCores)
	factorscores = foreach(omega = iter(omegamatrix, by='row', chunksize=chunksize),
						y = iter(Y, by='row', chunksize=chunksize),
						eta = iter(partialpredictor,by='row',chunksize=chunksize),
						.combine='rbind')	%dopar% {
		drawscoresBlock(y, loadings, omega, eta, overdisp,  qs, nuis,factorprecision)
	}
	return(factorscores)
}


drawscoresBlock = function(Y, loadings, omegamatrix, partialpredictor, overdisp, qs, nuis, factorprecision=NULL)
{
  N = nrow(Y)
	K = ncol(loadings)
  Eye = diag(1,K)
	if(missing(factorprecision)) factorprecision = diag(1,K)
	factorscores = matrix(0, nrow=N, ncol=K)
	for(i in 1:N)
	{
		kappa = (Y[i,] - overdisp)/2 - omegamatrix[i,]*partialpredictor[i,]
		postprec = crossprod(loadings, omegamatrix[i,]*loadings) + factorprecision
		L = t(chol(postprec))
		Li = forwardsolve(L, Eye)
		postvar = crossprod(Li)
		#postvar = chol2inv(U)
		#postvar = solve(crossprod(loadings, omega*loadings) + factorprecision)
		postmean = postvar %*% {crossprod(loadings, kappa)}
    #sparcity here
		prob1 = diag(dnorm(0,postmean,postvar))
    p1<-prob1 > 0.00001
		qhat<-rep(1,K)
    rhats<-dnorm(0,0,nuis)/prob1
		qhat[p1] <- rhats[p1]/((1-qs[p1])/qs[p1]+rhats[p1])
		
		factorscores[i,] = ifelse(runif(K)<qhat,postmean + crossprod(Li, rnorm(K)),0)
		#factorscores[i,] = rmvnorm(1, postmean, postvar)
	}
	return(factorscores)
}



### Main MCMC function
countfactormcmc = function(Y, K, overdisp=1, burn=50, nmc = 500, qs=NULL, verbose=NULL)
# Y is a (numsamples X numfeatures) matrix of count outcomes
# K is the number of factors
# overdisp is the first parameter in the NB parameterization, ie y ~ NB(overdisp, logit(linearfunc))
{
	if(missing(verbose)) verbose=nmc+burn
	
	# Bookkeeping
	P = ncol(Y)
	N = nrow(Y)
	numtrials = Y + overdisp
	
  nuis<-rep(0,K)
	# Set up placeholders to save the relevant draws
	save.loadings = array(0, dim=c(P,K,nmc))
	save.intercept = matrix(0, nmc, P)
	save.factormeans = matrix(0, nmc, N)
	save.factorscores = array(0, dim=c(N,K,nmc))
	
	# Initialize the chain
	loadings = matrix(rnorm(P*K, 0, 1), nrow=P, ncol=K)
	factorscores = matrix(rnorm(N*K,0,1), nrow=N, ncol=K)
	intercept = rep(0,P)
	factormeans = rep(0,N)
	
	# Main loop
	for(draw in 1:(nmc+burn))
	{
		if(draw %% verbose == 0) cat("\tDraw", draw, "\n")
		
		linearpredictor = cbind(intercept,loadings,1) %*% rbind(1,t(factorscores), factormeans)
		linearpredictor = t(linearpredictor)
		
		omega = drawomegamatrix(numtrials, linearpredictor)
		loadings = drawloadings(Y, factorscores, omega, factormeans, intercept, overdisp)
		intercept = drawintercept(Y, factorscores, omega, factormeans, loadings, overdisp)
		factormeans = drawfactormeans(Y, factorscores, omega, intercept, loadings, overdisp)
  
    ms<- apply(factorscores,2,function(x) {sum(abs(x)>.0000001)})
		#sjm get q for each column
		qs<- sapply( ms, function(m) {rbeta(1,1+m,1+N-m)})
		for(s in 1:K){
      nuis[s]<-  1/(rgamma(1,(c_s+ms[s])/2,(c_s*d_s+as.numeric(t(factorscores[,s])%*%factorscores[,s]))/2))      
		}
    
    factorscores = drawscores(Y, loadings, omega, factormeans, intercept, overdisp, qs, nuis)

		# save draws
		if(draw > burn)
		{
		  save.factorscores[,,draw-burn] = factorscores
			save.loadings[,,draw-burn] = loadings
			save.intercept[draw-burn,] = intercept
			save.factormeans[draw-burn,] = factormeans
		}
	}
	list(loadings = save.loadings, intercept = save.intercept, factormeans = save.factormeans, factorscores=save.factorscores)
}
