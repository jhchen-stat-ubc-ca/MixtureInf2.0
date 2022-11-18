#' rmix.exp
#'
#' @description This function generates the iid samples from the finite exponential mixture
#' @param n The sample size.
#' @param alpha A vector of the mixing probabilities.
#' @param theta A vector of the subpopulation means.
#'
#' @return
#' @export
#'
#' @examples
rmix.exp <- function (n, alpha, theta)  {
  m = length(alpha)
  alpha = alpha/sum(alpha)
  data=c()
  
  nindex=rmultinom(1, n, alpha)
  for( i in 1:m) data=c(data,rexp(nindex[i],1/theta[i]))
  sample(data)
}