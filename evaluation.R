#####Evaluation Script########
#load("result66.Rdata")
fulldata = read.csv("speachMemcovariateDTM.csv")  #2310 speaches by 8468 words + covariates (in columns c(2:12, 8480) )
doc66 = which(fulldata$year67 == 0)
words  = c(colnames(fulldata[,13:ncol(fulldata)-1]),"party2c","type.xcounty","type.university","miniTRUE","malaportionment")
words = words[2:length(words)]

PV.mean = read.csv("PVmean.csv",row.names = 1)  ## mean proportion of variance P*K matrix
lmB66 = read.csv("lmB66.csv",row.names = 1)     # posterior B matrix P*K
lmF66 = read.csv("lmF66.csv",row.names = 1)     # posterior F matrix N*K
mean.dimwords = read.csv("meandimwords.csv",row.names = 1) # mean dim for words P-vector

#
#BBn66 = result66$BBn;  # Store B matrix from Ndrawn     P*K*Ndrawn
#FFn66 = result66$FFn;  # Store F matrix from Ndrawn     N*K*Ndrawn
#Alfn66 = result66$Alfn; # Store alpha vector from Ndrawn   P*Ndrawn
#Gamn66 = result66$Gamn; # Store gamma vector from Ndrawn  N*Ndrawn
#P = dim(BBn66)[1]  # 8467 words + 4 logit + 1 cont
#K = dim(BBn66)[2]
#Ndrawn = dim(BBn66)[3]
#N = dim(FFn66)[1]

##  posterior loading matrix Bhat
if(0){
lmB66 = matrix(0,P,K)
for(i in 1:P){
  for(j in 1:K){
    lmB66[i,j] = mean(BBn66[i,j,])
  }
}

## posterior factor matrix Fhat
lmF66 = matrix(0,N,K)
for(i in 1:N){
  for(j in 1:K){
    lmF66[i,j] = mean(FFn66[i,j,])
  }
}

}
#### Produce P*P Correlation Matrix postcov66   
#### postcov66[i,j] show correlation between word i and word j
if(1){
postBBt66 = lmB66%*%t(lmB66)
diagD66 = 1/sqrt(diag(postBBt66))
postcov66 = postBBt66
for(j in 1:P){
    postcov66[j,] = (postBBt66[j,]*diagD66[j]*diagD66)
}
}

## show most correlated words to "leasehold"
cov1 = sort(postcov66[which(words == "leasehold"),],decreasing = TRUE, index.return = T)
MostCov1 = words[cov1$ix[1:20]]
MostCov1

library(calibrate)
#### 1. Where are the speeches?
plot(lmF66[,6],lmF66[,7],xlab = "F6", ylab = "F7", main = "where are the speeches")
textxy(lmF66[,6],lmF66[,7],doc66)
abline(h = 0, v = 0, col = "gray60")

#### 2. Where are the words?
plot(lmB66[,2],lmB66[,6],xlab = "F2", ylab = "F6", main = "where are the words")
textxy(lmB66[,2],lmB66[,6],words)
abline(h = 0, v = 0, col = "gray60")

#### 4. dim of words
if(0){
dimwords = matrix(0,P,Ndrawn)
for(i in 1:Ndrawn){
  for(j in 1:P){
    dimwords[j,i] = length(which(abs(BBn66[j,,i])>0.000000000001))
  }
}


mean.dimwords = apply(dimwords,1,mean)
}

plot(1:8467,mean.dimwords[1:8467])
textxy(1:8467,mean.dimwords[1:8467],words[1:8467])


#### 3. Highest loading words in Factor 2
cov1 = sort(lmB66[,2],decreasing = TRUE, index.return = T)
MostCov1 = words[cov1$ix[1:20]]
MostCov1

## 5. variational proportion
if(0){
PV = BBn66
for(ind in 1:Ndrawn){
  for(i in 1:P){
    temp = sum(BBn66[i,,ind]*BBn66[i,,ind])
    ### in case total variance is tooo small
    if(temp<0.000000001){
      PV[i,,ind] = rep(0,K)
    }
    else{
      for(j in 1:K){
        PV[i,j,ind] = BBn66[i,j,ind]^2/temp
        
      }
    }
  }
  
}

PV.mean = matrix(0,P,K)
for(i in 1:P){
  for(j in 1:K){
    PV.mean[i,j] = mean(PV[i,j,])
  }
}

}

##plot PV
factoruse = PV.mean[1:8467,2]
plot(1:8467,factoruse,ylab = "proportion of variation",xlab = "word",main = "Factor 6")
textxy(1:8467,factoruse,words[1:8467])








