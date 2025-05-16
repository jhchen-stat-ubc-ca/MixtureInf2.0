#' rmix.beta
#'
#' @description This function generates the iid samples from the finite beta mixture
#' @param n The sample size.
#' @param mix_prop A vector of the mixing proportions.
#' @param alpha A vector of the component alphas.
#' @param beta A vector of the component betas.
#'
#' @return It returns a sample of size n from an m-component beta mixture.
#'
#' @examples n = 1000
#' mix_prop = c(.5, .3, .2)
#' alpha = c(1,20,3)
#' beta = c(10,0.5,2)
#' rmix.beta(n, mix_prop, alpha,beta)
#' @export
rmix.beta <- function(n, mix_prop, alpha, beta) {
  nindex <- as.vector(rmultinom(1, n, mix_prop))
  data <- unlist(mapply(rbeta, n = nindex, shape1 = alpha, shape2 = beta, SIMPLIFY = FALSE))
  sample(data)
}
