rm(list=ls())


source("functions.R")
source("SCselection.R")
library(MASS)
library(parallel)
library(FBFsearch)
# generate true data
n <- 100
p <- 200
cat("n",n,"p",p,"\n")

# set sparsity level 4% or 8% (percentage of nonzero entries in true L)
spars <- 4
cat("spars",spars,"\n")

# first create a matrix with all entries being 1
true.L <- matrix(1,p,p)

# set the upper triangular to be 0
true.L[upper.tri(true.L)] <- 0

# generate and select nonzero entries in each column of L such that the expected pertenage of nonzero etries in each column of L is the pre-determined sparsity value
for(i in 1:(p-1))
{
  true.L[(i+1):p,i] <- true.L[(i+1):p,i]*rbinom(p-i,1,min(spars/100,1))
}
nonzero <- as.vector(1*colSums(true.L))
for(j in 1:(p-1)) {
  if(nonzero[j] == 1) true.L[p,j] = 0.5
}


# create the adjacency matrix according to the true L
true.adj <- true.L
diag(true.L) <- 1
diag(true.adj) <- 0

# generate the true covariance matrix and data based on the true covariance matrix
true.Sigma <- solve(true.L%*%t(true.L))
X <- mvrnorm(n,rep(0, p),Sigma=true.Sigma)
S <- 1/n*t(X)%*%X
tau <- 1
print(tau)


# set the initial value for the MCMC algorithm
tildeS <- S + diag(p)/(n*tau)
epsilon <- 0.3
ch <- chol(solve(S+epsilon*diag(p))) 
dd <- diag(ch)
L.epsilon <- t(ch/sqrt(dd))
L.epsilon1 <- L.epsilon
for(i in 2:p){
  for(j in 1:(i-1)){
    if(abs(L.epsilon[i,j]) < 0.1) L.epsilon1[i,j] <- 0
    else L.epsilon1[i,j] <- 1}}
nonzero.adj <- as.vector(1*colSums(L.epsilon1!=0))
for(m in 1:(p-1)) {
  if(nonzero.adj[m] == 1) L.epsilon1[p,m] = 1
}
diag(L.epsilon1) <- 0

# set hyperparameters
r <- 2
alpha1 <- 0.01
alpha2 <- 0.01
niter <- 10000
nburn <- 5000


##########non-local parallell
xuan = function(S, n, r, alpha1, alpha2, tau) {
  Zout <- matrix(0,p,p)
  p = nrow(S)
  Z <- L.epsilon1
  update_col = function(j) {
    # out = rep(0, p)
    
    
    
    curr_Zj = Z[(j+1):p, j]
    sum(curr_Zj)
    if(sum(curr_Zj) == 0) curr_Zj[sample(1:length(curr_Zj),1)] = 1
    curr_post = logposterior(curr_Zj,j)
    
    for (i in 1:niter) {
      cand_Zj = curr_Zj
      if (runif(1) < 0.5) {
        idx = sample(which(cand_Zj > 0), 1)
        cand_Zj[idx] = 0
      } else {
        idx = sample(which(cand_Zj == 0), 1)
        cand_Zj[idx] = 1
      }
      if(sum(cand_Zj) == 0) cand_Zj[sample(1:length(cand_Zj),1)] = 1
      cand_post = logposterior(cand_Zj,j)
      if (log(runif(1)) < min(cand_post-curr_post, 0)) {
        curr_Zj = cand_Zj
        curr_post = cand_post
      }
      if(i > nburn) Zout[,j] <- Zout[,j] + c(rep(0, j), curr_Zj)
      
      if (i %% 10000 == 0) {
        cat("column",j,i, "iterations completed.\n")
      }
    }
    
    return(Zout[,j])
  }
  
  mclapply(1:(p-10), update_col, mc.cores = 1)
}

########tau_1
time = Sys.time()
haha <- xuan(S, n, r, alpha1, alpha2, tau)
elapsed.nonlocal = Sys.time() - time
print(elapsed.nonlocal)

# evaluate the method performance
dag.nonlocal <- cbind(1*(matrix(unlist(haha), nrow = p)/(niter-nburn)>0.5),L.epsilon1[,((p-9):p)])
evaluation.dag(true.adj, dag.nonlocal)







#######tau_2
tau <- 2
tildeS <- S + diag(p)/(n*tau)

time = Sys.time()
haha1 <- xuan(S, n, r, alpha1, alpha2, tau)
elapsed.nonlocal = Sys.time() - time
print(elapsed.nonlocal)


dag.nonlocal <- cbind(1*(matrix(unlist(haha1), nrow = p)/(niter-nburn)>0.5),L.epsilon1[,((p-9):p)])
evaluation.dag(true.adj, dag.nonlocal)












#######tau_p
tau <- p^2.01
tildeS <- S + diag(p)/(n*tau)

time = Sys.time()
haha1 <- xuan(S, n, r, alpha1, alpha2, tau)
elapsed.nonlocal = Sys.time() - time
print(elapsed.nonlocal)


dag.nonlocal <- cbind(1*(matrix(unlist(haha1), nrow = p)/(niter-nburn)>0.5),L.epsilon1[,((p-9):p)])
evaluation.dag(true.adj, dag.nonlocal)






















##############FBF search
time = Sys.time()
M_q=FBF_LS(cor(X), n, matrix(0,p,p), 0, 0.01, 5000)
elapsed.FBF = Sys.time() - time
print(elapsed.FBF)
G_med=M_q
G_med[M_q>=0.5]=1
G_med[M_q<0.5]=0 #median probability DAG
evaluation.dag(true.adj,t(G_med))












#################Lasso with quantile tuning
lambdas <- c(1,2*qnorm(p=(0.1/(2*p*(1:(p-1)))),lower.tail=FALSE))/sqrt(n)
LassoAdjseq <- DAGLassoseq(Y=X,lambda.seq=lambdas)
LassoAdjseq[which(LassoAdjseq!=0)] <- 1
evaluation.dag(true.adj, LassoAdjseq)









################SC (Lee) search
# hyperparameters
alpha = 0.999
gamma = 1
nu0 = 0.1
c1 = 0.0001
c2 = 2
niter = 5000
nburn = 0.2*niter
nadap = 0



time = Sys.time()
res = SCselection(X, alpha, gamma, nu0, c1, c2, c3=NULL, init.A=NULL, niter, nburn, nadap, sample.dj=FALSE)
elapsed.lee = Sys.time() - time
print(elapsed.lee)


num.err = rep(0, p)
num.one = rep(0, p)
num.correct.one = rep(0, p)
num.false.one = rep(0, p)
p0.bar = 0
p1.bar = 0
incl.pr <- matrix(0,p,p)
for(j in 2:p){
  Sj.mat = as.matrix(res[[j]]$Sj.mat)
  incl.pr[j,seq(1:(j-1))] <- apply(Sj.mat, 2, mean)}
print(evaluation.dag(true.adj,1*(incl.pr>0.5)))











