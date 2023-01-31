#' Top one method for a vector
#'
#' This function allows you to identify the largest entry of a vector, set it to 1 and set other entries to 0.
#' @param x A vector
#' @return A vector with entries as 1 or 0.
#' @export
top_one = function(x){
  ind.max=which.max(x)
  x[ind.max]=1;x[-ind.max]=0
  return(x)}