#' Select the optimal threshold to sparsify the estimated orthogonal mapping matrix. 
#'
#' This function applies cross-validation to select the optimal threshold to sparsify the estimated orthogonal mapping matrix.
#'
#' @param ys a vector for outcome in source data (if any).
#' @param xs the design matrix for source data.
#' @param yt a vector for outcome in target data (if any).
#' @param xt the design matrix for target data.
#' @param family a description of the error distribution and link function to be used in the model. The default is gaussian.
#' @param nfold the number of folds in cross-validation
#' @param refine whether we refine the estimated mapping matrix leveraging the first moment relationship before sparsifing. The default is not to refine the estimated mapping matrix. 
#' @param threshold_vec a vector of candidate thresholds
#' @return The optimal threshold that can be used to sparsify the mapping matrix.
#' @export
cv_modtransfer = function(ys, xs, yt, xt, family = 'gaussian', nfold = 5, refine = F, 
                          threshold_vec = c(0.00001,0.0001,0.001,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
  partition_s <- createFolds(ys, k = nfold, list = TRUE, returnTrain = FALSE)
  partition_t <- createFolds(yt, k = nfold, list = TRUE, returnTrain = FALSE)
  p = ncol(xs)
  pred.trans.gen.all = matrix(NA,nrow = nfold, ncol = length(threshold_vec))
  for (t in 1:length(threshold_vec)) {
    threshold_tmp = threshold_vec[t]
    for (partid in 1:nfold) {
      ind_s = partition_s[[partid]]
      ind_t = partition_t[[partid]]
      
      xs_tr = xs[-ind_s,]
      xt_tr = xt[-ind_t,]
      ys_tr = ys[-ind_s]
      yt_tr = yt[-ind_t]
      xs_test = xs[ind_s,]
      xt_test = xt[ind_t,]
      ys_test = ys[ind_s]
      yt_test = yt[ind_t]
      
      co_occurs = cov(xs_tr)
      co_occurt = cov(xt_tr)
      
      us = svd(co_occurs)$u
      ut = svd(co_occurt)$u
      
      #general method
      s = matrix(0, nrow = p, ncol = p)
      for (l in 1:p) {
        s[l,l] = ifelse(ut[,l]%*%us[,l]>0, 1, -1)
      }
      pi_tilde_gen = (us%*%s)%*%t(ut)
      
      if (refine){
        mu_t = apply(xt, 2, mean)
        mu_s = apply(xs, 2, mean)
        pi_hat_gen = matrix(NA, nrow = p, ncol = p)
        for (m in 1:p) {
          pi_hat_gen[m,] = pi_tilde_gen[m,] + ((mu_s[m]-pi_tilde_gen[m,]%*%mu_t)/norm(mu_t,'2')^2)*mu_t
        }
      }else{
        pi_hat_gen = pi_tilde_gen
      }
      
      pi_hat_gen_sparsify = pi_hat_gen
      pi_hat_gen_sparsify[abs(pi_hat_gen_sparsify)<threshold_tmp]=0
      
      coef.naive = coef(glm(ys_tr~-1+xs_tr, family = family))
      coef.trans.gen=(matrix(coef.naive,nrow=1)%*%pi_hat_gen_sparsify)
      
      if(family=='gaussian'){
        pred.trans.gen.all[partid,t] = norm(yt_test-xt_test%*%matrix(coef.trans.gen,nrow = p),"2")
      }else{
        pred.trans.gen.all[partid,t] = auc(yt_test,inv_logit(xt_test%*%matrix(coef.trans.gen,nrow = p)))
      }
    }
  }
  if(family=='gaussian'){
    threshold.opt = threshold_vec[which.min(apply(pred.trans.gen.all, 2, mean))]
  }else{
    threshold.opt = threshold_vec[which.max(apply(pred.trans.gen.all, 2, mean))]
  }
  return(threshold.opt)
}