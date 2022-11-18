#' qmix.exp
#'
#' @description This function computes the quantile of the finite exponential mixture distribution
#' @param alpha A vector of the mixing probabilities.
#' @param theta A vector of the subpopulation means.
#' @param qq The level of the quantile to be computed.
#'
#' @return
#' @export
#'
#' @examples
qmix.exp <- function(alpha, theta, qq) {
  qmix.exp.sub <- function(x,alpha,theta,qq) {
    sum(alpha*pexp(x, rate=1/theta))-qq }
 
  if((qq<0)|(qq>1)) stop("quantile level out of range")
  uniroot(qmix.exp.sub, c(0, qexp(qq, rate=1/max(theta))), 
          alpha = alpha, theta = theta, qq=qq, tol=1e-10)$root 
}