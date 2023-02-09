#' EMtest.norm.iter
#'
#' @param xx The input data that can be either a vector or a matrix with the 1st column being the observed values
#'          and the 2nd column being the corresponding frequency.
#' @param nn 
#' @param m0 The order of the finite mixture model.
#' @param para 
#' @param sigma0 
#' @param beta.i 
#' @param pen.size 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
EMtest.norm.iter <- function(xx, nn, m0, para, sigma0, beta.i, pen.size, k) {
  alpha = para[1:(2*m0)]
  mu = para[(2*m0+1):(4*m0)]
  sigma = para[(4*m0+1):(6*m0)];  sig = sigma^.5
  
  for (i in 1:k)  {     ## K EM-iterations with special penalty structure
    pdf.sub = xx.sq = c()
    for (j in 1:(2*m0)) {	
      pdf.sub= cbind(pdf.sub,dnorm(xx,mu[j],sig[j])*alpha[j])
      xx.sq = cbind(xx.sq, (xx-mu[j])^2)  }
    pdf = apply(pdf.sub, 1, sum) + 1e-50
    ww = pdf.sub/pdf; ww.tot = apply(ww,2,sum)
    
    mu = apply(ww*xx,2,sum)/ww.tot
    xx.sq.weighted = apply(ww*xx.sq, 2, sum)
    sigma = (xx.sq.weighted + 2*pen.size[2]*sigma0)/(ww.tot + 2*pen.size[2])
    alpha=c()
    for (j in 1:m0) {
      temp = ww.tot[2*j-1]+ww.tot[2*j]
      if (ww.tot[2*j-1]/temp <=0.5)
        beta.i[j] = min((ww.tot[2*j-1]+pen.size[1])/(temp + pen.size[1]), 0.5)
      else
        beta.i[j] = max((ww.tot[2*j-1])/(temp + pen.size[1]), 0.5)
      alpha = c(alpha, c(beta.i[j], (1-beta.i[j]))*temp/nn) }
  }
  
  ###Compute the penalized log-likelihood value
  ln = sum(log(dmix.norm(xx, alpha, mu, sigma^.5)+ 1e-50))
  pen1 = pen.size[1]*sum(log(1-abs(1-2*beta.i)))
  temp = array(rbind(sigma0,sigma0))/sigma
  pen2 = -pen.size[2]*sum(temp - log(temp)-1)
  pln = ln + pen1 + pen2
  
  output = c(alpha, mu, sigma, pln)
}