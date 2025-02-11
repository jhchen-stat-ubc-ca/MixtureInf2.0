#' pmle.beta.sub
#'
#' @description A sub function for pmle.beta, does the actual work of PMLE of the beta mixture.
#' It is used in the pmle.beta function.
#' @export
pmle.beta.sub <- function(x,m0,para0,an,epsilon) {
  mix_porp = para0[1:m0]
  alpha = para0[(m0+1):(2*m0)]
  beta = para0[(2*m0+1):(3*m0)]
  theta = para0[(m0+1):(3*m0)]
  # E step
  n = length(x)
  pdf.sub = t(mapply(function(pp, aa, bb) pp * dbeta(x, aa, bb), mix_porp, alpha, beta))
  pdf.mixture = apply(pdf.sub,2,sum) +1e-100
  # M step
  ww = sweep(pdf.sub, 2, pdf.mixture, FUN = "/")
  mix_porp = (rowSums(ww)+epsilon)/(n+m0*epsilon)
  theta = Pen_M_Step(x, t(ww), mix_porp, theta, an)
  alpha = theta[1:m0]
  beta = theta[(m0+1):(2*m0)]
  loglike = sum(log(dmix.beta(x, mix_porp, alpha, beta)+1e-100))
  ploglik=sum(log(dmix.beta(x, mix_porp, alpha, beta))+1e-100)+
    sum(log(alpha*exp(-alpha))+log(beta*exp(-beta)))/an
  ind = sort(alpha,index.return = TRUE)$ix
  return(c(mix_porp[ind], alpha[ind],beta[ind], loglike,ploglik))
}


#' Pen_M_Step
#'
#' @description A sub function for pmle.beta.sub, does the actual work of the M step in the EM algorithm.
#' It is used in the pmle.beta.sub function.
#' @export
Pen_M_Step <- function(x, ww, mix_porp, theta, an) {
  ploglikelihood <- function(theta, x) {
    if (any(theta <= 0)) {return(NA)}
    n = length(x)
    alpha = theta[1:length(mix_porp)]
    beta = theta[(length(mix_porp)+1):(2*length(mix_porp))]
    return(ww*log(dmix.beta(x, mix_porp, alpha, beta))+(log(alpha*exp(-alpha))+log(beta*exp(-beta)))/an)
  }
  pgradlik <- function(theta,x) {
    if (any(theta <= 0)) {return(NA)}
    alpha = theta[1:length(mix_porp)]
    beta = theta[(length(mix_porp)+1):(2*length(mix_porp))]
    g1j <- sapply(1:length(alpha), function(j) {
      sum(ww[,j]*(log(x) + digamma(alpha[j] + beta[j]) - digamma(alpha[j])))-((alpha[j]-1)/alpha[j])/an
    })
    
    g2j <- sapply(1:length(beta), function(j) {
      sum(ww[,j]*(log(1 - x) + digamma(alpha[j] + beta[j]) - digamma(beta[j])))-((beta[j]-1)/beta[j])/an
    })
    return(c(g1j,g2j))
  }
  phesslik <- function(theta,x) {
    if (any(theta <= 0)) {return(NA)}
    alpha = theta[1:length(mix_porp)]
    beta = theta[(length(mix_porp)+1):(2*length(mix_porp))]
    k = length(alpha)
    n = length(x)
    Hess = matrix(0,ncol=2*length(alpha),nrow=2*length(alpha))
    for (j in 1:k) {
      Hess[j, j] <- sum(ww[,j]*(trigamma(alpha[j] + beta[j])- trigamma(alpha[j])))-(1/alpha[j]^2)/an
      Hess[j + k, j + k] <- sum(ww[,j]*(trigamma(alpha[j] + beta[j]) - trigamma(beta[j])))-(1/beta[j]^2)/an
      Hess[j, j + k] <- sum(ww[,j]*trigamma(alpha[j] + beta[j]))
      Hess[j + k, j] <- Hess[j, j + k]
    }
    return(Hess)
  }
  new_est_param = maxLik::maxNR(fn=ploglikelihood,grad=pgradlik,hess=phesslik,
                        start = theta, x = x)
  return(new_est_param$estimate)
}