load("01-1803CorpDTM4546.RData")

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

# number of documents
N = dim(YY)[1];

# number of terms
P = dim(YY)[2];

## number of topics 
Ntopic = 5;
K = Ntopic;

# prior

## Construct Prior Here
V = diag(1, K)
V0 = diag(1,K)
a = 0.01
niu_a = 10
niu_alpha = 10
niu = rep(0.98,K)
q = rep(0.9,K)
h = 1
c_s = 1
d_s = 1

testPrior = list(V,V0,a,niu_a,niu_alpha,niu,q,h,c_s,d_s);
names(testPrior) = c("V","V0","a","niu_a","niu_alpha","niu","q","h","c_s","d_s");

## Construct Initial for Gibbs Sampling  
## Can start with the final iteration of last run
## B = result$BBn[,,Ndraw]...

B = matrix(rnorm(P*K,0,1), nrow=P, ncol=K)
f = matrix(rnorm(N*K,0,1), nrow=N, ncol=K)
alpha = rnorm(P, 0, niu_alpha)
gamma = rnorm(N, 0, 1)

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




