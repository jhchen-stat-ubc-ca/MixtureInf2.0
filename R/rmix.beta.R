#' rmix.beta
#'
#' @description This function generates the iid samples from the finite beta mixture
#' @param n The sample size.
#' @param porp A vector of the mixing proportions.
#' @param alpha A vector of the component alphas.
#' @param beta A vector of the component betas.
#'
#' @return It returns a sample of size n from an m-component beta mixture.
#'
#' @examples n = 1000
#' porp = c(.5, .3, .2)
#' alpha = c(1,20,3)
#' beta = c(10,0.5,2)
#' rmix.beta(n, porp, alpha,beta)
#' @export
rmix.beta <- function (n, porp, alpha, beta)  {
  m = length(porp)
  data=c()
  
  nindex=rmultinom(1, n, porp)
  for( i in 1:m) data=c(data,rbeta(nindex[i],alpha[i],beta[i]))
  sample(data)
}
