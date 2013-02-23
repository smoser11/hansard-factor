### nohup Rscript reformIllustration.R  > NBf_6667reformParty.log 2>&1 &

#illustration of NBfactor model using floor speeches, HoC 1866-1867

setwd("~/Dropbox/research/nbfactor-hansard")
rm(list=ls())

load("~/Dropbox/research/Hansard rhetoric/memberID-plainText-DTM-date-party-speechLevel.RData")

#load("reformIllSpeechSpeaker.RData")
# speakerDTM66<-data.frame(speakerDTM66)
# speakerDTM67<-data.frame(speakerDTM67)
# speakerDTM66$year<-1866
# speakerDTM67$year<-1867
# speakerDTM66$memberID<-row.names(speakerDTM66)
# speakerDTM67$memberID<-row.names(speakerDTM67)
# speakerDTM<-merge(speakerDTM66,speakerDTM67,all=T,incomparables=NA)
# save(speakerDTM66, speakerDTM67, speakerDTM, speechLevel_allNospeech,file="reformIllSpeechSpeaker.RData")

DTM<-comb.dat[,-c(1:5,7)] #col 6 = 'cons' 0if liberal, 1 if conservative

library(MASS)
library(tm)
library(coda)
source("countfaMCMC-foreach-sjmEdit.R")

#  Parameter Setting
# How many topics?
K = 10
niu = rep(0.98,K)
q = rep(0.9,K)
h = 1
c_s = 2  #prior on nui
d_s = 1  #prior on nui , currently only enters in qs, sparcity
# Negative-binomial hyperparameter
overdisp = 1


# Data set properties
N = dim(DTM)[1]
P = dim(DTM)[2]
vocab<-colnames(DTM)  

# Load document-term matrix
Y = as.matrix(DTM)
  
## Run the MCMC and extract the summaries
mcmc1 = countfactormcmc(Y, K=K, nmc=2, burn=1, verbose=1)
  
any(mcmc1$factorscores[,,]==0)

bhat = apply(mcmc1$loadings,c(1,2), mean)
ahat = colMeans(mcmc1$intercept)
muhat = colMeans(mcmc1$factormeans) 

fileStem="6667reformK5"
eval(parse(text=paste("save(mcmc1,bhat,ahat,muhat, file=\"out-",fileStem,".RData\") ",sep="") ))
  
dir()[grep(dir(),pattern="RData")]
load("out-6667reformParty.RData")

length(vocab)
dim(mcmc1$loadings)


#check convervence -- 'convergence.R'


#visualize factors

bhat = apply(mcmc1$loadings,c(1,2), mean)
ahat = colMeans(mcmc1$intercept)
muhat = colMeans(mcmc1$factormeans)

for(m in 1:5000){
rownames(mcmc1$loadings[,,m])<-vocab
}

rownames(bhat)<-vocab  


#for each factor, k, get the top.words highest loading vocab words

top.words=15;

factors.table<-{}
require(xtable)
for(k in 1:K) {
 factors.table<-cbind(factors.table, vocab[order(bhat[,k],decreasing=TRUE)[1:top.words]] 
) 
}

#sort columbs by 'least cons' to 'most cons'
factors.table[,order(bhat["cons",])]

xtable(factors.table[,order(bhat["cons",])])

# descriptive statistics
##speeches, total, by party year

require(lattice)
#histogram of speech lengths
histogram(rowSums(comb.dat[,-c(1:7)]) #col 6 = 'cons' 0if liberal, 1 if conservative
)

pdf(file="histSpeechDist.pdf")
hist(rowSums(comb.dat[,-c(1:7)]),breaks=100 #col 6 = 'cons' 0if liberal, 1 if conservative
,xlab="Speech length (word count)", main="Histogram of Speech Length"  )
dev.off()

speechDesc<-data.frame(memberID=comb.dat$memberID, 
                       date=as.factor(ifelse(as.Date(comb.dat$V1)<=as.Date("26-jun-1866","%d-%b-%Y"),"1866","1867" )) ,
                       party=as.factor(ifelse(comb.dat$cons==1,"cons","lib")), length=rowSums(comb.dat[,-c(1:7)]))

summary(speechDesc)
xtable(summary(speechDesc))
var(speechDesc$length)

var(speechDesc$length[speechDesc$length!=max(speechDesc$length)])

table(speechDesc$date,speechDesc$party)


# look at changes over time:

##load model fit on 66, 67 data seperately
rm(list=ls())

load("out-66reformYearParty.RData")
load("out-67reformYearParty.RData")
rownames(bhat66)<-vocab
rownames(bhat67)<-vocab  

#10 most popular words in 66,67
vocab[order(ahat66,decreasing=TRUE)[1:11]]
vocab[order(ahat67,decreasing=TRUE)[1:11]]


##average 'dimension' of a speech: in 66 and in 67
dim(mcmc1$loadings)
any(mcmc1$factorscores[,,]==0)


##change in `discrimination' of words in 66 to 67 - variance of 


##partisanship over time

bhat66["cons",]
bhat67["cons",]

summary(bhat66["cons",])
var(bhat66["cons",])

summary(bhat67["cons",])
var(bhat67["cons",])

pdf(file="partyAffect.pdf")
par(mfrow=c(1,2))
scatter.smooth(main="1866",xlab="Variance, party loadings",ylab="Mean, party loadings",  span=1/3,pch=".",degree=2,apply(mcmc66$loadings[1,,],2, var ), apply(mcmc66$loadings[1,,],2, mean ) )
scatter.smooth(main="1867",xlab="Variance, party loadings", ylab="",ylim=c(-.3,.3),  span=1/3,pch=".",degree=2,apply(mcmc67$loadings[1,,],2, var ), apply(mcmc67$loadings[1,,],2, mean ) )
dev.off()

scatter.smooth(xlab="Variance, party loadings",ylab="Mean, party loadings", main="Posterior Mean/Var, Conservative Party \n in 1866 (K=10)", span=1/3,pch=".",degree=2,apply(mcmc66$loadings[1,,],2, var ), apply(mcmc66$loadings[1,,],2, mean ) )
scatter.smooth(xlab="Variance, party loadings",ylab="Mean, party loadings", main="Posterior Mean/Var, Conservative Party \n in 1867 (K=10)", span=1/3,pch=".",degree=2,apply(mcmc67$loadings[1,,],2, var ), apply(mcmc67$loadings[1,,],2, mean ) )

     
#density plot of varaiance of party loading in 66, 67
plot(density(apply(mcmc66$loadings[1,,],2, var )),main="Discrimination of party,\n legislative discourse 1866",xlab="Variance of party loadings, 10 factors")
plot(density(apply(mcmc67$loadings[1,,],2, var )),main="Discrimination of party,\n legislative discourse 1867",xlab="Variance of party loadings, 10 factors")


w=vocab[grep(vocab,pattern="suffra")]
w
which(vocab==w)
w %in% vocab
scatter.smooth(xlab="Variance, word loadings",ylab="Mean word loading", main="Posterior Mean/Var, ``suffra'' \n in 1866 (K=10)", span=1/3,pch=".",degree=2,apply(mcmc66$loadings[which(vocab==w),,],2, var ), apply(mcmc66$loadings[which(vocab==w),,],2, mean ) )
scatter.smooth(xlab="Variance, word loadings",ylab="Mean word loading", main="Posterior Mean/Var, ``suffra'' \n in 1867 (K=10)", span=1/3,pch=".",degree=2,apply(mcmc67$loadings[which(vocab==w),,],2, var ), apply(mcmc67$loadings[which(vocab==w),,],2, mean ) )


vocab[grep(vocab,pattern="constitu")]
w=vocab[grep(vocab,pattern="constit")][1]
w
which(vocab==w)
w %in% vocab
pdf(file="constitu.pdf")
par(mfrow=c(1,2))
scatter.smooth(xlab="Variance, word loadings",ylab="Mean word loading", main=" ``constitu'' \n in 1866", span=1/3,pch=".",degree=2,apply(mcmc66$loadings[which(vocab==w),,],2, var ), apply(mcmc66$loadings[which(vocab==w),,],2, mean ) )
scatter.smooth(xlab="Variance, word loadings",ylab="Mean word loading", main=" ``constitu'' \n in 1867", span=1/3,pch=".",degree=2,apply(mcmc67$loadings[which(vocab==w),,],2, var ), apply(mcmc67$loadings[which(vocab==w),,],2, mean ) )
dev.off()

vocab[grep(vocab,pattern="princip")]
w=vocab[grep(vocab,pattern="princip")][2]
w
which(vocab==w)
w %in% vocab
pdf(file="princip.pdf")
par(mfrow=c(1,2))
scatter.smooth(xlab="Variance, word loadings",ylab="Mean word loading", main=" ``princip'' \n in 1866", span=1/3,pch=".",degree=2,apply(mcmc66$loadings[which(vocab==w),,],2, var ), apply(mcmc66$loadings[which(vocab==w),,],2, mean ) )
scatter.smooth(xlab="Variance, word loadings",ylab="Mean word loading", main=" ``princip'' \n in 1867", span=1/3,pch=".",degree=2,apply(mcmc67$loadings[which(vocab==w),,],2, var ), apply(mcmc67$loadings[which(vocab==w),,],2, mean ) )
dev.off()


scatter.smooth(xlab="Variance, word loadings",ylab="Mean, word loadings", main="Posterior Mean/Var, ``princip'' \n in 1866 (K=10)", span=1/3,pch=".",degree=2,apply(mcmc66$loadings[which(vocab==w),,],2, var ), apply(mcmc66$loadings[which(vocab==w),,],2, mean ) )
scatter.smooth(xlab="Variance, word loadings",ylab="Mean, word loadings", main="Posterior Mean/Var, ``princip'' \n in 1867 (K=10)", span=1/3,pch=".",degree=2,apply(mcmc67$loadings[which(vocab==w),,],2, var ), apply(mcmc67$loadings[which(vocab==w),,],2, mean ) )




plot(density(apply(mcmc66$loadings[which(vocab==w),,],2, var )),main="Discrimination  1866",xlab="Variance of factor score")
plot(density(apply(mcmc67$loadings[which(vocab==w)
,,],2, var )),main="Discrimination  1867",xlab="Variance of factor score")



##compare with LDA-like topic models

ldaReformDTM<-dtm2ldaformat(as.simple_triplet_matrix(as.matrix(DTM)), omit_empty = TRUE)

library(tm)
library(topicmodels)
library(lda)

result <- lda.collapsed.gibbs.sampler(ldaReformDTM$documents,
  K,  ## Num clusters
  ldaReformDTM$vocab,
  500,  ## Num iterations
  0.1,
  0.1,
  compute.log.likelihood=TRUE) 

top.words=10

top.words.v <- top.topic.words(result$topics, top.words, by.score=TRUE)

require(xtable)
xtable(top.words.v)



top.topic.words(result$topics,num.words=15)

topic.proportions <- t(result$document_sums) / colSums(result$document_sums)

sample.docs<- c(347, 283) #gladstone resignation of ministry, col. lindsay representaino of people bill vote    sample(1:dim(topic.proportions)[1], 6)

topic.proportions <-
  topic.proportions[sample.docs,]

topic.proportions[is.na(topic.proportions)] <-  1 / K

colnames(topic.proportions) <- apply(top.words.v, 2, paste, collapse=" ")

topic.proportions.df <- melt(cbind(data.frame(topic.proportions),
                                   document=factor(1:length(sample.docs))),
                             variable.name="topic",
                             id.vars = "document")  

qplot(topic, value, fill=document, ylab="proportion",
      data=topic.proportions.df, geom="bar") +
        opts(axis.text.x = theme_text(angle=90, hjust=1)) +  
        coord_flip() +
        facet_wrap(~ document, ncol=dim( topic.proportions.df)[1]/K)








s.result <- slda.em(documents=ldaReformDTM$documents,
K=10,      vocab=ldaReformDTM$vocab,
             num.e.iterations=100,
             num.m.iterations=40,
             alpha=1.0, eta=.1, params=params,
                    annotations=comb.dat$cons,
           variance=0.25,
           lambda=1.0,
             logistic=FALSE,
           method="sLDA")

require("ggplot2")

 Topics <- apply(top.topic.words(s.result$topics, 5, by.score=TRUE),
                  +                 2, paste, collapse=" ")

 coefs <- data.frame(coef(summary(s.result$model)))
theme_set(theme_bw())

 coefs <- cbind(coefs, Topics=factor(Topics, Topics[order(coefs$Estimate)]))

 coefs <- coefs[order(coefs$Estimate),]

 qplot(Topics, Estimate, colour=Estimate, size=abs(t.value), data=coefs) +
   geom_errorbar(width=0.5, aes(ymin=Estimate-Std..Error,
                   ymax=Estimate+Std..Error)) + coord_flip()








