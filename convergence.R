#convergence diagnostics of NB sparse Factor

library(coda)

load( ~/Dropbox/research/nbfactor-hansard/out-1811CorpDTM5634.RData)

load(~/Dropbox/research/NBfactor/out-6667reformK5.RData)

factormeans.mc<-mcmc(mcmc1$factormeans)
intercept.mc<-mcmc(mcmc1$intercept)


#reshpae factor loadings into NUM_ITTER by K*V matrix
library(R.utils)
library(matlab)

dim(mcmc1$loadings)
loadingsLong<-matlab::reshape(mcmc1$loadings, dim(mcmc1$loadings)[1]*dim(mcmc1$loadings)[2],dim(mcmc1$loadings)[3]   )
loadingsLong<-t(loadingsLong)  #now in NUM_itter by VocabSize * NumFactors

loadings.mc<-mcmc(loadingsLong)

plot(intercept.mc[,sample(1:dim(intercept.mc)[2],10)])
plot(loadings.mc[,sample(1:dim(loadings.mc)[2],10)])
plot(factormeans.mc[,sample(1:dim(factormeans.mc)[2],10)])


heidel.diag(loadings.mc)
raftery.diag(loadings.mc[,10])
raftery.diag(intercept.mc[,100])

geweke.diag(intercept.mc)
geweke.diag(loadings.mc)
#how many factor-loadings (cells) pass the gweke difference in means test?
countGeweke.diag <- function(mc) {
c( sum( sapply( apply(mc,2,geweke.diag)  ,function(x) { x$z<1.96 })),dim(mc)[2])
}
countGeweke.diag(intercept.mc[,sample(1:dim(intercept.mc)[2],1000 )])
countGeweke.diag(loadings.mc[,sample(1:dim(loadings.mc)[2],1000 )  ])

countHeidel.diag <- function(mc) {
  a<- apply(mc,2,heidel.diag)
  c(sum(abs(a[6,]/a[5,]) <.1 ,na.rm=TRUE),length(na.omit(a[6,]/a[5,])),dim(mc)[2])
}

countHeidel.diag(intercept.mc[,sample(1:dim(intercept.mc)[2],100)])
countHeidel.diag(loadings.mc[,sample(1:dim(loadings.mc)[2],100)])

heidel.diag(loadings.mc[,sample(1:dim(loadings.mc)[2],100)])
