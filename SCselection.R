# X : n X p data matrix
SCselection <- function(X, alpha, gamma, nu0, c1, c2, c3=NULL, init.A=NULL, niter, nburn, nadap, sample.dj=FALSE){
   
   n = nrow(X)
   p = ncol(X)
   if(!exists("lars")) library(lars)
   res = list()
   
   for(j in 2:p){
      
      Sj.mat = matrix(0, nrow=niter + nburn, ncol=j-1)
      sj.mat = rep(0, niter + nburn)
      Rsq.mat = rep(0, niter + nburn)
      logpost.mat = rep(0, niter + nburn)
      
      # set Rj value
      if(is.null(c3)){
         Rj = floor( n/(log(p, base=10)*log(n, base=10)) )
      }else{
         Rj = floor( n/log(p)*min(c3, 1/log(n)) )
      }
      
      # data (tilde.Xj) and design matrix (Z.j) for j-th dimension
      tilde.Xj = as.matrix(X[, j])
      Z.j = as.matrix(X[, 1:(j-1)])
      tilde.Xj.prod = sum( tilde.Xj^2 )
      
      
      # initialzation
      if(is.null(init.A)){
         # initial guess for Sj based on lasso
         ob.j = lars(x = Z.j, y = tilde.Xj, normalize = FALSE, intercept = FALSE, use.Gram = FALSE)
         cv.j = cv.lars(x = Z.j, y = tilde.Xj, plot.it = FALSE, se = FALSE, normalize = FALSE, intercept = FALSE, use.Gram = FALSE)
         las.est = coef(ob.j, s = cv.j$index[which.min(cv.j$cv)], mode = "fraction")
         if(sum(las.est != 0) > min(j-1, Rj)){
            las.est[ sample(which(las.est != 0), size = sum(las.est != 0) - min(j-1, Rj), replace = FALSE) ] = 0
         }
         Sj.init = as.numeric(las.est != 0)
      }else{
         Sj.init = as.numeric(init.A[j, 1:(j-1)] != 0)
      }
      
      if(sum(Sj.init) == 0) Sj.init[sample(1:(j-1), size = 1)] = 1
      
      
      # initial guess for 'dj' (if it is needed)
      if(sample.dj){ # if we also want to infer dj's
         dj.mat = rep(0, niter + nburn)
         # res.prod.beta = as.numeric( tilde.Xj - Z.j[, Sj.init > 0] %*% las.est[Sj.init > 0] )
         # dj.init = sum(res.prod.beta^2)/max(n-sum(Sj.init), 1) # initial guess for dj (based on lasso)
         dj.init = 1/rgamma(n = 1, shape = nu0/2, rate = 10) # initial guess for dj (based on prior)
         dj = dj.init
         dj.mat[1] = dj
      }else{
         dj = NULL
      }
      
      
      Sj = Sj.init
      logpost.old = logpost(X, j, Rj, Sj, dj, alpha, gamma, nu0, c1, c2)
      
      if(is.finite(logpost.old$val) == FALSE){
         # initial guess for Sj based on lasso
         ob.j = lars(x = Z.j, y = tilde.Xj, normalize = FALSE, intercept = FALSE, use.Gram = FALSE)
         cv.j = cv.lars(x = Z.j, y = tilde.Xj, plot.it = FALSE, se = FALSE, normalize = FALSE, intercept = FALSE, use.Gram = FALSE)
         las.est = coef(ob.j, s = cv.j$index[which.min(cv.j$cv)], mode = "fraction")
         if(sum(las.est != 0) > min(j-1, Rj)){
            las.est[ sample(which(las.est != 0), size = sum(las.est != 0) - min(j-1, Rj), replace = FALSE) ] = 0
         }
         Sj.init = as.numeric(las.est != 0)
         Sj = Sj.init
         logpost.old = logpost(X, j, Rj, Sj, dj, alpha, gamma, nu0, c1, c2)
      }
      
      Sj.mat[1, ] = Sj
      sj.mat[1] = sum(Sj)
      Rsq.mat[1] = 1 - n*logpost.old$dj / tilde.Xj.prod
      logpost.mat[1] = logpost.old$val
      logpost.uniq.mat = logpost.old$val
      dj.uniq.mat = logpost.old$dj
      unique.index = 1
      logpost.new = list()
      
      accept.ind = 0 # number of acceptance in MH sampler
      for(i in 2:(niter+nburn)){
         
         if(i <= nadap){
            Sj.new = Sprop.adap(X, Sj, j, Rj) # new proposal of Sj (Sj.new)
            log.prop.kernel = lkernelprob.adap(X, Sj, Sj.new, j) # log proportion of proposal kernel
         }else{
            Sj.new = Sprop(Sj, j, Rj) # new proposal of Sj (Sj.new)
            log.prop.kernel = lkernelprob(Sj, Sj.new, j) # log proportion of proposal kernel
         }
         
         # calcultaion of logpost.new
         if(sample.dj){ # if we also want to infer dj's
            logpost.new = logpost(X, j, Rj, Sj.new, dj, alpha, gamma, nu0, c1, c2)
            match.ind = 1
         }else{ # no inference for dj's
            
            if(is.na(unique.index[2])){
               match.ind = prod(as.numeric(Sj.mat[1,] == Sj.new))
            }else{
               match.ind = compare.to.rows(as.matrix(Sj.mat[unique.index ,]), Sj.new)
            }
            
            if(match.ind == 0){
               logpost.new = logpost(X, j, Rj, Sj.new, dj, alpha, gamma, nu0, c1, c2)
            }else{
               logpost.new$val = logpost.uniq.mat[match.ind]
               logpost.new$dj = dj.uniq.mat[match.ind]
            }
            
         }
         
         # MH sampler for Sj index
         if(is.finite(logpost.new$val) == FALSE) logpost.new$val = -Inf
         if(runif(1) <= exp(logpost.new$val - logpost.old$val + log.prop.kernel)){
            Sj = Sj.new
            logpost.old = logpost.new
            if(match.ind == 0){
               logpost.uniq.mat = c(logpost.uniq.mat, logpost.old$val)
               dj.uniq.mat = c(dj.uniq.mat, logpost.old$dj)
               unique.index = c(unique.index, i)
            }
            accept.ind = accept.ind + 1
         }
         
         # sampling for 'dj' (if it is needed)
         if(sample.dj){
            dj = 1/rgamma(1, shape=(alpha*n+nu0)/2, rate=alpha*n*dShat(X, j, Sj)/2)            
            dj.mat[i] = dj
         }
         
         Sj.mat[i, ] = Sj
         sj.mat[i] = sum(Sj)
         Rsq.mat[i] = 1 - n*logpost.old$dj / tilde.Xj.prod
         logpost.mat[i] = logpost.old$val
         
      } # end of (i in 2:(niter+nburn)) for loop
      
      res[[j]] = list(Sj.mat = Sj.mat[-(1:nburn),], sj.mat = sj.mat[-(1:nburn)], Rsq.mat = Rsq.mat[-(1:nburn)], logpost.mat = logpost.mat[-(1:nburn)])
     # cat("Posterior sampling for ", j,"th row is completed. . . . . .")
    #  cat(accept.ind, " samples are accepted.\n")
      
   } # end of (j in 2:p) for loop
   
   return(res)
}


###################################################################
# Auxiliary functions
###################################################################

compare.to.rows <- function(SS, S){
   h <- function(v) sum(abs(v - S))
   o <- apply(SS, 1, h)
   if(all(o > 0)){
      return(0)
   }else{
      return(which(o == 0)[1])
   }
}


Sprop <- function(S, j, Rj){
   s = sum(S)
   upper.ind = min(j-1, Rj)
   
   if(s == 0){ # if current S has no index
      S[sample(which(S == 0), 1)] = 1
   }else if(s == upper.ind){ # if current S has maximum index (upper.ind)
      S[sample(which(S > 0), 1)] = 0
   }else{
      
      if(runif(1) <= 0.5){ # introducing additional one 0 to current S
         if(s == 1){
            S[which(S == 1)] = 0
         }else{
            S[sample(which(S > 0), 1)] = 0
         }
         
      }else{ # introducing additional one 1 to current S
         if(j-1-s == 1){
            S[which(S == 0)] = 1
         }else{
            S[sample(which(S == 0), 1)] = 1
         }
      }
      
   }
   return(S)
}
Sprop.adap <- function(X, S, j, Rj){
   n = nrow(X)
   s = sum(S)
   upper.ind = min(j-1, Rj)
   
   if(s == 0){ # if current S has no index
      S[sample(which(S == 0), 1)] = 1
   }else if(s == upper.ind){ # if current S has maximum index (upper.ind)
      S[sample(which(S > 0), 1)] = 0
   }else{
      
      if(runif(1) <= 0.5){ # introducing additional one 0 to current S
         if(s == 1){
            S[which(S == 1)] = 0
         }else{
            prob.of.new.ind = rep(0, s)
            for(k in 1:s){
               S.old = S
               S.old[which(S == 1)[k]] = 0
               prob.of.new.ind[k] = 1/dShat(X, j, S.old)
            }
            S[sample(which(S == 1), 1, prob = prob.of.new.ind)] = 0
         }
      }else{ # introducing additional one 1 to current S
         if(j-1-s == 1){
            S[which(S == 0)] = 1
         }else{
            prob.of.new.ind = rep(0, j-1-s)
            for(k in 1:(j-1-s)){
               S.old = S
               S.old[which(S == 0)[k]] = 1
               prob.of.new.ind[k] = 1/dShat(X, j, S.old)
            }
            S[sample(which(S == 0), 1, prob = prob.of.new.ind)] = 1
         }
      }
      
   }
   return(S)
}


lkernelprob <- function(Sj, Sj.new, j){
   if(sum(Sj - Sj.new) == -1){
      return(log(j-1 - sum(Sj)) - log(sum(Sj.new)))
   }else{
      return(log(sum(Sj)) - log(j-1 - sum(Sj.new)))
   }
}
lkernelprob.adap <- function(X, Sj, Sj.new, j){
   if(sum(Sj - Sj.new) == -1){
      prob.of.new.ind = rep(0, sum(Sj==0))
      for(k in 1:sum(Sj==0)){
         S.old = Sj
         S.old[which(Sj == 0)[k]] = 1
         prob.of.new.ind[k] = 1/dShat(X, j, S.old)
      }
      q.denom = prob.of.new.ind[which(Sj == 0) == which(Sj - Sj.new == -1)] / sum(prob.of.new.ind)
      prob.of.new.ind = rep(0, sum(Sj.new==1))
      for(k in 1:sum(Sj.new==1)){
         S.old = Sj.new
         S.old[which(Sj.new == 1)[k]] = 0
         prob.of.new.ind[k] = 1/dShat(X, j, S.old)
      }
      q.numer = prob.of.new.ind[which(Sj.new == 1) == which(Sj.new - Sj == 1)] / sum(prob.of.new.ind)
      
      return(log(q.numer) - log(q.denom))
   }else{
      prob.of.new.ind = rep(0, sum(Sj==1))
      for(k in 1:sum(Sj==1)){
         S.old = Sj
         S.old[which(Sj == 1)[k]] = 0
         prob.of.new.ind[k] = 1/dShat(X, j, S.old)
      }
      q.denom = prob.of.new.ind[which(Sj == 1) == which(Sj - Sj.new == 1)] / sum(prob.of.new.ind)
      prob.of.new.ind = rep(0, sum(Sj.new==0))
      for(k in 1:sum(Sj.new==0)){
         S.old = Sj.new
         S.old[which(Sj.new == 0)[k]] = 1
         prob.of.new.ind[k] = 1/dShat(X, j, S.old)
      }
      q.numer = prob.of.new.ind[which(Sj.new == 0) == which(Sj.new - Sj == -1)] / sum(prob.of.new.ind)
      
      return(log(q.numer) - log(q.denom))
   }
}


logpost <- function(X, j, Rj, S, dj, alpha, gamma, nu0, c1, c2){
   # If dj=NULL, it returns the log of marginal posterior prob. \pi(S | data)
   # If dj!=NULL, it returns the log of conditional posterior prob. \pi(S | dj, data)
   n = nrow(X)
   p = ncol(X)
   s = sum(S)
   
   dS.hat = dShat(X, j, S)
   
#    logpij = -s*log(c1) - c2*s*log(p) - lchoose(j-1, s) + log(s <= min(j-1, Rj))
   logpij = dpois(s, lambda = 1, log = T) - lchoose(j-1, s) + log(s <= min(j-1, Rj)) # Poisson prior for the model size
   
   if(is.null(dj)){
   	  logpost.val = logpij - s*log(1 + alpha/gamma)/2 - (alpha*n + nu0)*log(dS.hat)/2
      return( list(val = logpost.val, dj = dS.hat) )
   }else{
      logpost.val = logpij - s*log(1 + alpha/gamma)/2 - alpha*n*dS.hat/(2*dj)
      return( list(val = logpost.val, dj = dS.hat) )
   }
}


dShat <- function(X, j, S){
   n = nrow(X)
   s = sum(S)
   tilde.Xj = as.matrix(X[, j])
   
   if(s == 0){
      return( sum((tilde.Xj)^2)/n )
   }else{
      Zj = as.matrix(X[, 1:(j-1)])
      X.Sj = as.matrix(Zj[, S>0])
      VhatX = sum(tilde.Xj^2)/n   # scalar
      Covhat = t(X.Sj)%*%matrix(tilde.Xj)/n   # s X 1 matrix
      res = VhatX - t(Covhat)%*%solve( t(X.Sj)%*%X.Sj/n )%*%Covhat
      return( res )
   }
}

