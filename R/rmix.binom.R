#' rmix.binom
#'
#' @description  This function generates a random sample from the finite binomial mixture distributions.
#' @param n The sample size.
#' @param size The number of trials.
#' @param alpha A vector of the mixing proportions.
#' @param theta A vector of the probabilities of success of each component.
#'
#' @return
#' @export
#'
#' @examples #Generate a random sample from a 3 component binomial mixture,
#' and computes the sample mean and variance. 
#' alpha = c(.5, .1, .4)
#' theta = c(.1, .35, .8)
#' size = 25
#' x = rmix.binom(1000, size, alpha, theta)
#' mean(x)
#' var(x)
rmix.binom <- function(n, size, alpha, theta) {
  if(any(alpha<0)) stop("negative mixing probabilities")
  if(any(theta<0)) stop("negative subpopulation probabilities")
  m = length(alpha)
  alpha =alpha/sum(alpha)
  data =c()
  nindex=rmultinom(1,n,alpha)
  for(i in 1:m)
    data=c(data,rbinom(nindex[i],size,theta[i]))
  sample(data)   
  ## avoid separation of the observations from different subs.
}