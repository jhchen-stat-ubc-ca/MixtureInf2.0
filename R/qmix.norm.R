#' qmix.norm
#'
#' @description This function computes the quantile of the univariate finite normal mixture one by one.
#' It is used in the tildeB22.norm0 function.
qmix.norm <- function(pp, alpha, mu, sigma, tol=1e-8) {
  if(any(alpha<0))
    stop("error: negative mixing proportion")
  if(any(sigma<0))
    stop("error: negative standard deviation")
  if(length(pp) > 1)
    stop("compute one quantile at a time")
  if(pp < 0| pp> 1)
    stop("pp out of (0, 1)")
  m1 = length(alpha)
  m2 = length(mu)
  m3 = length(sigma)
  if((m1-m2)^2+(m1-m3)^2 > 0)
    stop("error: differ lengths of alpha, mu and sigma")
  alpha = alpha/sum(alpha)
  tt = rep(0, m1)
  for(i in 1:m1) tt[i] = qnorm(pp, mu[i], sigma[i])
  low = min(tt);  high = max(tt)
  n.iter = (log(1/tol)*(high-low))
  for(i in 1:n.iter) {
    tmp = (low+high)/2
    if(pmix.norm(tmp, alpha, mu, sigma) > pp) {
      high = tmp} else {low = tmp}
  }
  return((low+high)/2)
}