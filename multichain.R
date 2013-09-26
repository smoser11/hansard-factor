library(CombFactor)
library(foreach)
library(iterators)
library(doMC)
registerDoMC()

load('data67setting.RData')

result67 = CombFactor(ContY67,Logit67,CountY67,testPrior,testInitial,10,200,100)

# result67 = foreach(i=1:8) %dopar% {
		# CombFactor(ContY67,Logit67,CountY67,testPrior,testInitial,10,2000,500)
	# }


tx = 1:ncol(result67$Alfn)
alphastart = apply(result67$Alfn, 1, function(x) lm(x ~ tx)$fitted.values[length(tx)])
