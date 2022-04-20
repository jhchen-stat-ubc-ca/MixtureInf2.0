#' qmixnorm
#'
#' @description This function computes the alp0-quantile of Normal mixture.
#' @param alpha A vector of the mixture probabilities. 
#' @param theta A vector of means from each component.
#' @param alp0 A given probability.	
#'
#' @return
#' @export
#'
#' @examples
qmixnorm <-
  function(alpha,theta,alp0)
  {
    uniroot(pmixnorm,c(qnorm(alp0,min(theta)),qnorm(alp0,max(theta))),alpha=alpha,theta=theta,alp0=alp0)$root 
  }