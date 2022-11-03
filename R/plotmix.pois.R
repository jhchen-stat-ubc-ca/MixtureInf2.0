## x: a vector of data whose histogram is to be drawn.
## x.max: if x is null, 0:x.max is the range of pmf to be plotted.
##  alpha: a vector of mixing proportions
##  mu : a vector of subpopulation means
## sub.pmf: logical, draw subpopulation pmf's.
## nclass: number of bars for histogram.

#' plotmix.pois
#'
#' @description This function plots the data in histogram with 
#'              the corresponding fitted mixture pmf
#'              or the mixture pmf only, for Poisson mixture. 
#'              The pmf is given as if it is a continuous curve.
#' @param x A vector of data whose histogram is to be drawn from.
#' @param x.max If x is null, 0:x.max is the range of pmf to be plotted.
#' @param alpha A vector of the mixing proportions.
#' @param theta A vector of the subpopulation means.
#' @param sub.pmf logical, draw subpopulation pmf's.
#' @param nclass The number of bars for the histogram.
#' @param main The title for the histogram.
#' @param xlab The name of the x-axis label.
#' @param ylab The name of the y-axis label.
#'
#' @return
#' @export
#'
#' @examples n = 1000
#' mu = c(9, 10)
#' alpha = c(.4, .6)
#' x = rmix.pois(200, alpha, mu)
#' plotmix.pois(x, alpha=alpha, mu = mu)
plotmix.pois = function(x = NULL, x.max= NULL, alpha, mu, 
                        sub.pmf = T, nclass=NULL, main="", 
                        xlab="Counts", ylab="Prob") {
  if(is.null(x)) {
    pmf.x = dmix.pois(0:x.max, alpha, mu)
    plot(0:x.max, pmf.x, type="l", main=main, xlab = xlab, ylab = ylab)
    for(i in 1:length(alpha)) {
      pmf.x = dpois(0:x.max, mu[i])*alpha[i]
      lines(0:x.max, pmf.x, lty = 2) }
  }
  
  if(is.matrix(x)) {
    xx = c()
    for(i in 1:dim(x)[1]) xx = c(xx, rep(x[i,1], x[i,2]))
    x = as.numeric(xx)
  }
  
  if(is.vector(x)) {
    hist(x, nclass = nclass, probability = T,
         main=main, xlab = xlab, ylab = ylab)
    pmf.x = dmix.pois(1:max(x), alpha, mu)
    lines(x = pmf.x, col = "red")
    if(sub.pmf) {
      for(i in 1:length(alpha)) {
        pmf.x = dpois(1:max(x), mu[i])*alpha[i]
        lines(x = pmf.x, col = "blue")  } 
    } 
  }
}