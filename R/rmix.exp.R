#' rmix.exp
#'
#' @description This function generates the iid samples from the finite exponential mixture
#' @param n The sample size.
#' @param alpha A vector of the mixing probabilities.
#' @param theta A vector of the subpopulation means.
#'
#' @return It returns a sample of size n from an m-component exponential mixture.
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#'
#' @examples n = 4000
#' mu = c(3, 9, 18)
#' alpha = c(.5, .3, .2)
#' rmix.exp(n, alpha, mu)
#' @export
rmix.exp <- function (n, alpha, theta)  {
  m = length(alpha)
  alpha = alpha/sum(alpha)
  data=c()
  
  nindex=rmultinom(1, n, alpha)
  for( i in 1:m) data=c(data,rexp(nindex[i],1/theta[i]))
  sample(data)
}