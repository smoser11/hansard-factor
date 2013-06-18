# load data
load("01-1803CorpDTM4546.RData")

# load package
library(NBfactor)

## preprocessing YY: doc*term matrix
DTM = reformDTM
NR = as.numeric(DTM[4])
NC = as.numeric(DTM[5])
totl = length(DTM[1]$i)
YY = matrix(0,NR,NC)
for(index in 1:totl){
  YYi = as.numeric(DTM[1]$i[index])
  YYj = as.numeric(DTM[2]$j[index])
  YY[YYi,YYj] = as.numeric(DTM[3]$v[index])
}
## YY(i,j) is the number of occurances of word j in document i 

# number of documents
N = dim(YY)[1];

# number of terms
P = dim(YY)[2];

## number of topics 
Ntopic = 5;
K = Ntopic;

# Model:
# YY_ij  ~ NB(h, exp(Psi_ij)/1+exp(Psi_ij)) i = 1...N, j = 1...p
# Psi_i ~ N(alpha + gamma_i*1 +    B*f_i, sigma^-1) i = 1...N
#  p*1     p*1             p*1  p*k k*1    p*p

## Construct Prior Here
V = diag(1, K)        # f_i ~ N(0, V) prior for factor
V0 = diag(1,K)   
a = 0.01              # gamma_i ~ N(a, 1)
niu_a = 10            # a~ N(0, niu_a)
niu_alpha = 10        # alpha_j ~ N(0, niu_alpha)
niu = rep(0.98,K)     # niu_s~ IG(c_s, c_s*d_s/2)
q = rep(0.9,K)        # b_js ~ q_s*N(0,niu_s) + (1-q_s)*delta0   sparsity
h = 1
c_s = 1
d_s = 1

testPrior = list(V,V0,a,niu_a,niu_alpha,niu,q,h,c_s,d_s);
names(testPrior) = c("V","V0","a","niu_a","niu_alpha","niu","q","h","c_s","d_s");

## Construct Initial for Gibbs Sampling  
## Can start with the final iteration of last run
## B = result$BBn[,,Ndraw]...

B = matrix(rnorm(P*K,0,1), nrow=P, ncol=K)    # loading matrix
f = matrix(rnorm(N*K,0,1), nrow=N, ncol=K)	  # factor matrix
alpha = rnorm(P, 0, niu_alpha)                # term-vise intercept
gamma = rnorm(N, 0, 1)                        # document-vise intercept

testInitial = list(B,f,alpha,gamma);
names(testInitial) = c("B","f","alpha","gamma");


## number of iterations to be burn
Nburn = 0;
## number of iterations to be drawn
Ndrawn = 1;


test1 = NBfactor(YY, testPrior,testInitial,Ntopic, Ndrawn, Nburn);

BBn = test1$BBn;  # Store B matrix from Ndrawn     P*K*Ndrawn
FFn = test1$FFn;  # Store F matrix from Ndrawn     N*K*Ndrawn
Alfn = test1$Alfn; # Store alpha vector from Ndrawn   P*Ndrawn
Gamn = test1$Gamn; # Store gamma vector from Ndrawn  N*Ndrawn



## evaluation

# plot paramter
if(0){
  #plot Gammar[par1]
  par1 = 4
  par(mfrow=c(3,1))
  ts.plot(Gamn[par1,],xlab="iterations",ylab="") 
  acf(Gamn[par1,],main="") # correlation
  hist(Gamn[par1,],prob=T,main="",xlab="")
}

if(0){
 #plot Alfar[par1]
  par1 = 15
  par(mfrow=c(3,1))
  ts.plot(Alfn[par1,],xlab="iterations",ylab="")
  acf(Alfn[par1,],main="")
  hist(Alfn[par1,],prob=T,main="",xlab="")
}

if(0){

 #calculate mean of drawn B's
  lmB = matrix(0,P,K)
  #index.drawn = (Nburn+1):Ndrawn
  #index.drawn=400:1000
  for(i in 1:P){
    for(j in 1:K){
      lmB[i,j] = mean(BBn[i,j,])
    }
  }

 # show the most popular terms in par2 th topic
  par2 = 1
  test1 = sort(lmB[,par2],decreasing=TRUE,index.return=T)

  MostTerm1=DTM[6][1]$dimnames$Terms[test1$ix[1:25]]
  print(MostTerm1)

}

# show alf: term-wise intercept
if(0){
  #index.drawn = 60:100  
  meanAlf = rep(0,P)
  for(i in 1:P){
  meanAlf[i] = mean(Alfn[i,])
  }
  test2 = sort(meanAlf,decreasing=TRUE,index.return=T)
  MostTerm2=DTM[6][1]$dimnames$Terms[test2$ix[1:25]]
  MostTerm2

}





