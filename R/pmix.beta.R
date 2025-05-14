#' pmix.beta
#'
#' @param x The cdf value of the beta mixture at x.
#' @param mix_porp A vector of the mixing proportions.
#' @param alpha A vector of the component alphas.
#' @param beta A vector of the component betas.
#'
#'
#' @examples 
#' x = c(rbeta(50,1,10),rbeta(50,9,0.5))
#' pmix_beta(x,c(.5,.5),c(1,9),c(10,0.5))
#' @export
pmix_beta <- function(x, mix_porp, alpha, beta) {
  cdf_matrix <- sapply(seq_along(mix_porp), function(i) {
    mix_porp[i] * pbeta(x, alpha[i], beta[i])
  })
  if (length(x) == 1) {
    sum(cdf_matrix)
  } else {
    rowSums(cdf_matrix)
  }
}
