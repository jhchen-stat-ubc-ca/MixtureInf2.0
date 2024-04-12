#' plotmix.binom
#'
#' @description This function plots a histogram together with a suggested mixture pmf.
#' @param x The input data, either a vector or a matrix with two columns: count and freq. 
#' @param size The number of trials.
#' @param alpha A vector of the mixing probabilities.
#' @param theta vector of probabilities of success of each component.
#' @param sub.pmf logical, draw subpopulation pmf's.
#' @param nclass The number of classes for a histogram.
#' @param main The main title of the output graph.
#' @param xlab The label name for the x-axis.
#' @param ylab The label name for the y-axis.
#'
#' @return It returns the histogram of the observations and the plot of the fitted probability mass function.
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#'
#' @examples alpha = c(.5, .1, .4)
#' theta = c(.1, .35, .8)
#' size = 25
#' x = rmix.binom(1000, size, alpha, theta)
#' plotmix.binom(x, size, alpha, theta, sub.pmf = T, nclass=size, main="", xlab="Counts", ylab="Prob")
#' @export
plotmix.binom <- function(x, size, alpha, theta, 
                          sub.pmf = T, nclass=NULL, main="", xlab="Counts", ylab="Prob") {
  if(is.null(x)) {
    pmf.x = dmix.binom(0:size, size, alpha, theta)
    plot(0:size, pmf.x, type="l", main=main, xlab = xlab, ylab = ylab)
    for(i in 1:length(alpha)) {
      pmf.x = dbinom(0:size, size, theta[i])*alpha[i]
      lines(0:size, pmf.x, lty = 2) }
  }
  
  if(is.matrix(x)) {
    xx = c()
    for(i in 1:dim(x)[1]) xx = c(xx, rep(x[i,1], x[i,2]))
    x = as.numeric(xx)
  }
  
  if(is.vector(x)) {
    pmf.x = dmix.binom(0:size, size, alpha, theta)
    hist(x+2, nclass = nclass, probability = T,
         main=main, xlab = xlab, ylab = ylab, ylim=c(0, max(pmf.x)*1.1))
    #hist(x+0.5, nclass = nclass, probability = T,
    #main=main, xlab = xlab, ylab = ylab, ylim=c(0, max(pmf.x)*1.1))
    lines(x = pmf.x, col = "red")
    if(sub.pmf) {
      for(i in 1:length(alpha)) {
        pmf.x = dbinom(0:size, size, theta[i])*alpha[i]
        lines(x = pmf.x, col = "blue")  } 
    } 
  }
}