fulldata = read.csv("speachMemcovariateDTM.csv")  #2310 speaches by 8468 words + covariates (in columns c(2:12, 8480) )
data66 = fulldata[which(fulldata$year67 == 0),]   # year 1866 data
CountY66 = as.matrix(data66[,13:(ncol(fulldata)-1)] )        # Count Data for NB factor
Logit66 = as.matrix(data66[,c(4,5,6,11)])                    # Logit Data for Logistic Factor
ContY66 = as.matrix(data66[,10])                             # Continuous Data For Factor


N = dim(CountY66)[1]
P = dim(CountY66)[2]
Pl = dim(Logit66)[2]
Pc = dim(ContY66)[2]

K = 10     # number of topic

## Construct Prior Here
V = diag(1, K)        # f_i ~ N(0, V) prior for factor
V0 = diag(1,K)   
a = 0              # gamma_i ~ N(a, 1)
niu_a = 10            # a~ N(0, niu_a)
niu_alpha = 10        # alpha_j ~ N(mean_alpha(j), niu_alpha)
mean_alpha = rep(0,P)
niu = rep(0.98,K)     # niu_s~ IG(c_s, c_s*d_s/2)
q = rep(0.98,K)        # b_js ~ q_s*N(0,niu_s) + (1-q_s)*delta0   sparsity
h = 1
c_s = 2
d_s = 1

a_l = 0;
niu_a_l = 10;
niu_alpha_l = 10;
niu_l = rep(10,K);

a_c = 0;
niu_a_c = 10;
niu_alpha_c = 10;
niu_c = rep(10,K);

nu0 = 3.0;
s02 = 1.0;



## Construct Initial for Gibbs Sampling  
## Can start with the final iteration of last run
## B = result$BBn[,,Ndraw]...

B = matrix(rnorm(P*K,0,1), nrow=P, ncol=K)    # loading matrix
f = matrix(rnorm(N*K,0,1), nrow=N, ncol=K)    # factor matrix
alpha = rnorm(P, 0, 1)                # term-vise intercept
gamma = rnorm(N, 0, 1)                        # document-vise intercept

B_l = matrix(rnorm(Pl*K,0,1), nrow=Pl, ncol=K)    # loading matrix
alpha_l = rnorm(Pl, 0, 1)                # term-vise intercept
gamma_l = rnorm(N, 0, 1)                        # document-vise intercept

B_c = matrix(rnorm(Pc*K,0,1), nrow=Pc, ncol=K)    # loading matrix
alpha_c = rnorm(Pc, 0, 1)                # term-vise intercept
gamma_c = rnorm(N, 0, 1)                        # document-vise intercept

testPrior = list(mean_alpha,a_c,niu_a_c,niu_alpha_c,niu_c,nu0,s02,a_l,niu_a_l,niu_alpha_l,niu_l,V,V0,a,niu_a,niu_alpha,niu,q,h,c_s,d_s);
names(testPrior) = c("mean_alpha","a_c","niu_a_c","niu_alpha_c","niu_c","nu0","s02","a_l","niu_a_l","niu_alpha_l","niu_l","V","V0","a","niu_a","niu_alpha","niu","q","h","c_s","d_s");

testInitial = list(B,f,alpha,gamma,B_l,alpha_l,gamma_l,B_c,alpha_c,gamma_c);
names(testInitial) = c("B","f","alpha","gamma","B_l","alpha_l","gamma_l","B_c","alpha_c","gamma_c");

result66 = CombFactor(ContY66,Logit66,CountY66, testPrior,testInitial,10,10000,0);