library(MASS)
#library(glmnet)

GetDiagSub <- function(ind,S1){
  return(S1[ind,ind])
}
GetVecSub <- function(ind,vec){
  return(vec[ind])
}
  SpecialDet <- function(x){
    if(length(x)==0) return(1)
    if(length(x)==1) return(as.numeric(x))
    return(det(x))
  }

evaluation.dag <- function(Adj1, Adj2){
  true.index <- which(Adj1==1)
  false.index <- which(Adj1==0)
  positive.index <- which(Adj2==1)
  negative.index <- which(Adj2==0)
  
  TP <- length(intersect(true.index,positive.index))
  FP <- length(intersect(false.index,positive.index))
  FN <- length(intersect(true.index,negative.index))
  TN <- length(intersect(false.index,negative.index))
  
  MCC.denom <- sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
  if(MCC.denom==0) MCC.denom <- 1
  MCC <- (TP*TN-FP*FN)/MCC.denom
  if((TN+FP)==0) MCC <- 1
  
  Precision <- TP/(TP+FP)
  if((TP+FP)==0) Precision <- 1
  FDR <- FP/(TP+FP)
  if((TP+FP)==0) FDR <- 1
  Recall <- TP/(TP+FN)
  if((TP+FN)==0) Recall <- 1
  Sensitivity <- Recall
  Specific <- TN/(TN+FP)
  if((TN+FP)==0) Specific <- 1

  return(list(FDR=FDR,Sensitivity=Sensitivity,MCC=MCC,FP=FP,TP=TP,FN=FN,TN=TN))
}


DAGLasso <- function(Y, lambda, maxitr=100, tol=1e-4){
  
  require(lassoshooting)  
  p = ncol(Y)
  n = nrow(Y)
  
  S = (t(Y) %*% Y)/n
  
  T = diag(p)
  D = rep(1, p)
  
  itr_log = eps_log = NULL
  
  for (k in 2:p){
    
    nuk_old = nuk_new = c(rep(0, k-1), 1)
    r = 0
    km1_ = 1:(k-1)
    
    repeat {
      
      r = r + 1
      
      
      nuk_new[k] = D[k]
      
      
      output = lassoshooting(XtX= S[km1_,km1_,drop=FALSE], 
                             Xty=-S[km1_,k], maxit = 100, lambda = 0.5*nuk_new[k]*lambda)
      nuk_new[km1_] = output$coefficients
      #print("hshsh")
      
      maxdiff = max(abs(nuk_new - nuk_old))
      if (maxdiff < tol || r >= maxitr){
        T[k, km1_] = nuk_new[km1_]
        eps_log = c(eps_log, maxdiff)
        itr_log = c(itr_log, r)
        break
      } else {
        nuk_old = nuk_new
      }
      
    }
    
  }
  return(T)
  
}



DAGLassoseq <- function(Y, lambda.seq, maxitr=100, tol=1e-4){
  require(lassoshooting)  
  p = ncol(Y)
  n = nrow(Y)
  #print("hahhaah")
  S = (t(Y) %*% Y)/n
  
  T = diag(p)
  D = rep(1, p)
  
  itr_log = eps_log = NULL
  
  for (k in 2:p){
    
    nuk_old = nuk_new = c(rep(0, k-1), 1)
    r = 0
    km1_ = 1:(k-1)
    
    repeat {
      
      r = r + 1
      
      
      nuk_new[k] = D[k]
      
      
      output = lassoshooting(XtX= S[km1_,km1_,drop=FALSE], 
                             Xty=-S[km1_,k], maxit = 100, lambda=0.5*nuk_new[k]*lambda.seq[k])
      nuk_new[km1_] = output$coefficients
      # print("hahhaah")
      
      maxdiff = max(abs(nuk_new - nuk_old))
      if (maxdiff < tol || r >= maxitr){
        T[k, km1_] = nuk_new[km1_]
        eps_log = c(eps_log, maxdiff)
        itr_log = c(itr_log, r)
        break
      } else {
        nuk_old = nuk_new
      }
      
    }
    
  }

  Adj <- matrix(0,p,p)
  Adj[which(T!=0)] <- 1
  for(i in 1:p) Adj[i,i] = 0
  return(Adj)
  
}


Evaluation <- function(beta1, beta2){
  true.index <- which(beta1==1)
  false.index <- which(beta1==0)
  negative.index <- which(beta2==0)
  
  TP <- length(intersect(true.index,positive.index))
  FP <- length(intersect(false.index,positive.index))
  FN <- length(intersect(true.index,negative.index))
  TN <- length(intersect(false.index,negative.index))
  
  
  Precision <- TP/(TP+FP)
  if((TP+FP)==0) Precision <- 1
  FDR <- FP/(TP+FP)
  if((TP+FP)==0) FDR <- 0
  Recall <- TP/(TP+FN)
  if((TP+FN)==0) Recall <- 1
  Sensitivity <- Recall
  Specific <- TN/(TN+FP)
  if((TN+FP)==0) Specific <- 1

  return(list(FDR=FDR,Sensitivity=Sensitivity,Specific=Specific,TP=TP,FP=FP,TN=TN,FN=FN))
}



log_posterior <- function(Zj, j){
  f = function(x, Zj) {
    LZj = x[1:sum(Zj)]
    dj = x[sum(Zj)+1]
    if(dj<0.5)
    {
      #cat("hahahahah")
      dj=10^5
    }
    #print(dj)
    Ij = c(rep(0, j), Zj)
    Sj = S[Ij, Ij]
    SZj = S[Ij, j]
    nz = sum(Zj)
    out = - nz*log(prod(seq(1, 2*r-1, 2))) - (nz+n)/2*log(2*pi) - (r*nz+(n+nz)/2+alpha1+1)*log(dj) -(r*nz+nz/2)*log(tau) -(n*(t(LZj)%*%Sj%*%LZj+S[j, j]+2*t(SZj)%*%LZj)/2/dj) - (t(LZj)%*%LZj/2/tau/dj+alpha2/dj) + 2*r*sum(log(abs(LZj)))
    return(-out)
  }
  
  posterior = function(obj, Zj) {
    nZj = sum(Zj)
    #print(nZj)
    Ij = c(rep(0, j), Zj)
    Sj = S[Ij, Ij]
    SZj = S[Ij, j]
    estimate = obj$estimate
    LZj_hat = estimate[1:nZj]
    dj_hat = estimate[nZj + 1]
    V11 = - diag(nZj) / tau / dj_hat - Sj *n / dj_hat - 2 * r / LZj_hat^2
    V12 = (Sj %*% LZj_hat + SZj) * n / dj_hat^2 + LZj_hat / tau / dj_hat^2
    V22 = (r * nZj + nZj / 2 + n / 2 + alpha1 + 1) / dj_hat^2 - sum(LZj_hat^2) / tau / dj_hat^3 - (t(LZj_hat) %*% Sj %*% LZj_hat + S[j, j] + 2 * t(SZj) %*% LZj_hat) * n / dj_hat^3 - 2 * alpha2 / dj_hat^3
    V = rbind(cbind(V11, V12), c(t(V12), V22))
    #print(-obj$minimum)
    return((nZj + 1) / 2 * log(2 * pi) - obj$minimum - sum(log(abs(eigen(V)$values)))/2)
  }
  
  estimates = nlm(f, c(L.epsilon[which(c(rep(0, j), Zj)!= 0),j], 1/dd[j]), Zj)
  post = posterior(estimates, Zj)
  return(post)
}



logposterior <- function(Zj,j) {
  Ij = c(rep(0, j), Zj)
  
  pa.set <- which(Ij!=0)
  pa.count <- length(pa.set)
  #print(pa.count)
  pa.mat.det <- 0
  pa.con <- 0
  logz <- 0
  if(pa.count == 0) logz = 0
  else{
    pa.mat.det = SpecialDet(GetDiagSub(pa.set, tildeS))
    pa.con = tildeS[j,j] - as.numeric(1*t(GetVecSub(as.vector(pa.set), tildeS[,j]))%*%solve(GetDiagSub(as.vector(pa.set), tildeS))%*%GetVecSub(pa.set, tildeS[,j]))
    logz = (-n/2-alpha1+pa.count/2)*log(pa.con*n/2 + alpha2) - 1/2*log(pa.mat.det) - (r+0.5)*pa.count*(log(tau)+log(n))+log(gamma((n-pa.count)/2+alpha1))
  }
  return(logz)
}
