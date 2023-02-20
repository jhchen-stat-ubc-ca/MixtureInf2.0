#' dmix.binom
#'
#' @description This function computes the probability mass function points for vector x.
#' It is used in various functions within the Binomial folder.
dmix.binom <- function(x, size, alpha, theta) {
  pmf.x = x*0
  m0 = length(alpha)
  for(i in 1:m0) { pmf.x = pmf.x + dbinom(x, size, theta[i])*alpha[i]}
  pmf.x
}
