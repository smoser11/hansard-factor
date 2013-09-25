library(CombFactor)
library(foreach)
library(iterators)
library(doMC)
registerDoMC()

load('data67setting.RData')


result67 = foreach(i=1:8) %dopar% {
		CombFactor(ContY67,Logit67,CountY67,testPrior,testInitial,10,2000,500)
	}
