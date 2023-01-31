#' Estimate the orthogonal mapping matrix that determines the heterogeneity between two heterogeneous datasets. 
#'
#' This function estimates the orthogonal mapping matrix leveraging the second moment relationship between the two datasets
#'
#' @param ys a vector for outcome in source data (if any).
#' @param xs the design matrix for source data.
#' @param yt a vector for outcome in target data (if any).
#' @param xt the design matrix for target data.
#' @param family a description of the error distribution and link function to be used in the model. The default is gaussian.
#' @param rowmax whether apply the top one method for each row or each column
#' @param PERMUATION whether one-to-one mapping or general mapping. The default is one-to-one mapping.
#' @param refine whether we refine the estimated mapping matrix leveraging the first moment relationship. The default is not to refine the estimated mapping matrix. 
#' @param sparsify whether we sparsify the estimated mapping matrix for general mappings with cross-validation. This step requires (partially observed) target outcomes.
#' @return A orthogonal mapping matrix.
#' @examples
#' n=30000
#' p=100
#' r=0.5
#' mu_predefine = rep(0,p)
#' A = matrix(runif(p^2)*2-1, ncol=p) 
#' cov.mat.predefine = t(A) %*% A
#' true.beta=rnorm(p)
#' x_s = MASS::mvrnorm(n = n, mu=mu_predefine, Sigma=cov.mat.predefine)
#' x_t = MASS::mvrnorm(n = n, mu=mu_predefine, Sigma=cov.mat.predefine)
#' 
#' k=round(p*r)
#' perm=c(sample(x=1:k,size=k,replace=F),(k+1):p)
#' perm.mat = matrix(rep(0,p*p),p)
#' for(j in 1:p){perm.mat[j,perm[j]]=1}
#' 
#' x_t_obs=x_t%*%perm.mat
#' y_s=x_s%*%true.beta+rnorm(n)
#' y_t=x_t%*%true.beta+rnorm(n)
#' 
#' rslt = modtran(ys=y_s, xs=x_s, yt=y_t, xt=x_t_obs, refine = F, sparsify = F)
#' 
#' coef.naive = coef(glm(y_s~-1+x_s, family = 'gaussian'))
#' coef.trans.gen=t(matrix(coef.naive,nrow=1)%*%rslt[["pi_hat_gen"]])
#' coef.trans.perm=t(matrix(coef.naive,nrow=1)%*%rslt[["pi_hat_perm"]])
#' 
#' norm(y_t-x_t_obs%*%coef.naive,"2")
#' norm(y_t-x_t_obs%*%coef.trans.perm,"2")
#' norm(y_t-x_t_obs%*%coef.trans.gen,"2")
#' @export
modtran = function(ys = NULL, xs, yt = NULL, xt, family = 'gaussian', rowmax = T, PERMUATION = T, refine = F, sparsify = F){
  if (ncol(xs) != ncol(xt)) 
    stop("Both systems must have the same number of features")
  
  if (PERMUATION & sparsify) 
    warning("We will always sparsify estimated mapping matrix under permutation-only mappings")
  
  #derive embeddings via SVD
  p = ncol(xs)
  co_occurs = cov(xs)
  co_occurt = cov(xt)
  
  us = svd(co_occurs)$u
  ut = svd(co_occurt)$u
  
  #general method
  s = matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    s[i,i] = ifelse(ut[,i]%*%us[,i]>0, 1, -1)
  }
  pi_tilde_gen = (us%*%s)%*%t(ut)
  
  #permutation-only method 
  if (PERMUATION) {
    s = matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
      s[i,i] = ifelse(norm(sort(ut[,i])-sort(-us[,i]),'2') > norm(sort(ut[,i])-sort(us[,i]),'2'), 1, -1)
    }
    pi_tilde_perm = (us%*%s)%*%t(ut)
  }
  
  #refine pi_tilde_perm and pi_tilde_gen
  if (refine){
    mu_t = apply(xt, 2, mean)
    mu_s = apply(xs, 2, mean)
    pi_hat_gen = pi_hat_perm = matrix(NA, nrow = p, ncol = p)
    for (k in 1:p) {
      pi_hat_gen[k,] = pi_tilde_gen[k,] + ((mu_s[k]-pi_tilde_gen[k,]%*%mu_t)/norm(mu_t,'2')^2)*mu_t
      if (PERMUATION) {
        pi_hat_perm[k,] = pi_tilde_perm[k,] + ((mu_s[k]-pi_tilde_perm[k,]%*%mu_t)/norm(mu_t,'2')^2)*mu_t
      }
    }
  }else{
    pi_hat_gen = pi_tilde_gen
    if (PERMUATION) {
      pi_hat_perm = pi_tilde_perm
    }
  }
  
  if(PERMUATION){
    if(rowmax){#row max
      pi_hat_gen=t(apply(pi_hat_gen,1,top_one))
      pi_hat_perm=t(apply(pi_hat_perm,1,top_one))
    }else{#column max
      pi_hat_gen=t(apply(pi_hat_gen,2,top_one))
      pi_hat_perm=t(apply(pi_hat_perm,2,top_one))
    }
  }else{
    if(sparsify){
      opt.threshold = cv_modtran(ys, xs, yt, xt, family = family, refine = refine, nfold = 5, 
                                           threshold_vec = c(0.00001,0.0001,0.001,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
      pi_hat_gen[abs(pi_hat_gen) < opt.threshold] = 0
    }else{
      pi_hat_gen=pi_hat_gen
    }
  }
  
  ###return mapping matrix
  
  if(PERMUATION){
    rslt = list("pi_hat_gen"=pi_hat_gen,
                "pi_hat_perm"=pi_hat_perm)
  }else{
    rslt = list("pi_hat_gen"=pi_hat_gen)
  }
  
  
  return(rslt)
}