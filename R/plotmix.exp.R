#' plotmix.exp
#'
#' @description This function plots a histogram plus the corresponding fitted mixture pdf or the mixture pdf only, for a exponential mixture.
#' @param alpha A vector of the mixing probabilities.
#' @param mu A vector of the subpopulation means.
#' @param qq The range of the plot, which is the qq th quantile.
#' @param x The data whose histogram is to be drawn from, can either be a vector or a matrix.
#' @param logObs logObs=T if data is in logarithm transformation form.
#' @param sub.pdf logical, draw subpopulation pdf's if True.
#' @param h A scale value used by ylim.
#' @param nclass The number of bars for the histogram if specified.
#' @param main The title for the histogram.
#' @param xlab The name of the x-axis label.
#' @param ylab The name of the y-axis label.
#' @param ylim The range of the y-axis, used only for the histogram if specified.
#'
#' @return It returns the histogram of observations and the plot of the fitted density
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#'
#' @examples n = 4000
#' mu = c(3, 9, 18)
#' alpha = c(.5, .3, .2)
#' x = rmix.exp(n, alpha, mu) 
#' plotmix.exp(alpha, mu, qq = 0.995, x) 
plotmix.exp <- function(alpha, mu, qq = 0.995, x = NULL, logObs=F,
                        sub.pdf = T, h=1.1, nclass = NULL, main="", 
                        xlab = "Observations scaled", ylab = "Density", ylim=NULL) 
{
  if(!logObs) {
    if(is.null(x)) {
      x.max = qmix.exp(qq, alpha, mu)
      x.rang = seq(0, x.max, 0.005)
      pdf.x = dmix.exp(x.rang, alpha, mu)
      plot(x.rang, pdf.x, type="l", ylim=c(0, max(pdf.x)*h),
           main = main, xlab = NULL, ylab = ylab, col = "red")
      for(i in 1:length(alpha)) {
        pdf.x = dexp(x.rang, rate=1/mu[i])*alpha[i]
        lines(x.rang, pdf.x, lty = 2, col = "blue") }
    }
    
    if(is.matrix(x)) {
      xx = c()
      for(i in 1:dim(x)[1]) xx = c(xx, rep(x[i,1], x[i,2]))
      x = as.numeric(xx)
    }
    
    if(is.vector(x)) {
      hist(x, freq= F, nclass = nclass, main=main, xlab = xlab, ylab = ylab)
      x.rang = seq(0, max(x), 0.005)
      pdf.x = dmix.exp(x.rang, alpha, mu)
      lines(x.rang, pdf.x, col = "red")
      if(sub.pdf) {
        for(i in 1:length(alpha)) {
          pdf.x = dexp(x.rang, rate=1/mu[i])*alpha[i]
          lines(x.rang, pdf.x, col = "blue")  } 
      } 
    }
  }
  
  if(logObs) {
    if(is.null(x)) {
      x.max = qmix.exp(qq, alpha, mu)
      x.min = qmix.exp(1-qq, alpha, mu)
      x.rang = log(seq(x.min, x.max, 0.005))
      pdf.x = dmix.exp(x.rang, alpha, mu, logObs)
      plot(x.rang, pdf.x, type="l", ylim=c(0, max(pdf.x)*h),
           main = main, xlab = NULL, ylab = ylab, col = "red")
      for(i in 1:length(alpha)) {
        pdf.x = exp(x.rang-exp(x.rang)/mu[i])*alpha[i]/mu[i]
        lines(x.rang, pdf.x, lty = 2, col = "blue") }
    }
    
    if(is.matrix(x)) {
      xx = c()
      for(i in 1:dim(x)[1]) xx = c(xx, rep(x[i,1], x[i,2]))
      x = as.numeric(xx)
    }
    
    if(is.vector(x)) {
      hist(x, freq= F, nclass = nclass, main=main, xlab = xlab, ylab = ylab)
      x.rang = seq(min(x), max(x), 0.005)
      pdf.x = dmix.exp(x.rang, alpha, mu, logObs=T)
      lines(x.rang, pdf.x, col = "red")
      if(sub.pdf) {
        for(i in 1:length(alpha)) {
          pdf.x = exp(x.rang-exp(x.rang)/mu[i])*alpha[i]/mu[i]
          lines(x.rang, pdf.x, col = "blue")  } 
      } 
    }
  }
  
}