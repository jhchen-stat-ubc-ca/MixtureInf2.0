#' pmix.exp
#'
#' @description This function computes the cdf of the finite exponential mixture distribution
#' @param x The value at which we compute the cdf.
#' @param alpha A vector of the mixing probabilities.
#' @param theta A vector of the subpopulation means.
#'
#' @return
#' @export
#'
#' @examples
pmix.exp <- function(x, alpha, theta) {
  sum(alpha*pexp(x,rate=1/theta))
}