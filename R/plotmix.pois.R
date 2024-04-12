#' plotmix.pois
#'
#' @description This function plots the data in a histogram with 
#'              the corresponding fitted mixture pmf
#'              or the mixture pmf only, for Poisson mixture. 
#'              The pmf is given as if it is a continuous curve.
#' @param x A vector of data whose histogram is to be drawn from.
#' @param x.max If x is null, 0:x.max is the range of pmf to be plotted.
#' @param alpha A vector of the mixing proportions.
#' @param mu A vector of the subpopulation means.
#' @param sub.pmf logical, draw subpopulation pmf's.
#' @param nclass The number of bars for the histogram.
#' @param main The title for the histogram.
#' @param xlab The name of the x-axis label.
#' @param ylab The name of the y-axis label.
#'
#' @return It returns the histogram of the observations and the plot of the fitted probability mass function.
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#'
#' @examples n = 2000
#' mu = c(2, 9, 11)
#' alpha = c(.2, .3, .5)
#' xx = rmix.pois(n, alpha, mu)
#' plotmix.pois(alpha,mu,xx,extra.height = 1.05)
#' @export
plotmix.pois = function(alpha, mu, x = NULL, x.max= NULL,
                        sub.pmf = T, nclass=NULL, extra.height = 1.05, main="", 
                        xlab="Counts", ylab="Prob") {
  if(is.null(x)) {
    pmf.x = dmix.pois(0:x.max, alpha, mu)
    plot(0:x.max, pmf.x, type="l", main=main, xlab = xlab, ylab = ylab, col="red")
    for(i in 1:length(alpha)) {
      pmf.x = dpois(0:x.max, mu[i])*alpha[i]
      lines(0:x.max, pmf.x, lty = 2, col="blue") }
  }
  
  if(is.matrix(x)) {
    xx = c()
    for(i in 1:dim(x)[1]) xx = c(xx, rep(x[i,1], x[i,2]))
    x = as.numeric(xx)
  }
  
  if(is.vector(x)) {
    hist(x, nclass = nclass, probability = T, ylim = c(0,max((hist(x, freq = F, nclass = nclass))$density)*extra.height),
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