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

## number of iterations to be burn
Nburn = 0;
## number of iterations to be drawn
Ndrawn = 1;
## number of topics 
Ntopic = 5;

test1 = NBfactor(YY, Ntopic, Ndrawn, Nburn);

BBn = test1$BBn;  # Store B matrix from Ndrawn     P*K*Ndrawn
FFn = test1$FFn;  # Store F matrix from Ndrawn     N*K*Ndrawn
Alfn = test1$Alfn; # Store alpha vector from Ndrawn   P*Ndrawn
Gamn = test1$Gamn; # Store gamma vector from Ndrawn  N*Ndrawn




