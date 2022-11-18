#' plotmix.exp
#'
#' @description This function plots a histogram plus the corresponding fitted mixture pdf or the mixture pdf only, for exponential mixture.
#' @param alpha A vector of the mixing probabilities.
#' @param theta A vector of the subpopulation means.
#' @param qq The range of the plot, which is the qq th quantile.
#' @param x The data whose histogram is to be drawn, can either be a vector or a matrix.
#' @param sub.pdf logical, draw subpopulation pdf's if True.
#' @param h A scale value used by ylim.
#' @param nclass The number of bars for the histogram if specified.
#' @param ss 
#' @param main The title for the histogram.
#' @param xlab The name of the x-axis label.
#' @param ylab The name of the y-axis label.
#' @param ylim Used only for the histogram if specified.
#'
#' @return
#' @export
#'
#' @examples
plotmix.exp <- function(alpha, theta, qq = 0.995, x = NULL, 
                        sub.pdf = T, h=1.1, nclass = NULL, ss = 30, main="", 
                        xlab = "Observations scaled", ylab = "Density", ylim=NULL) 
{
  if(is.null(x)) {
    x.max = qmix.exp(alpha, theta, qq)
    x.rang = seq(0, x.max, 0.005)
    pdf.x = dmix.exp(x.rang, alpha, theta)
    plot(x.rang, pdf.x, type="l", ylim=c(0, max(pdf.x)*h),
         main = main, xlab = NULL, ylab = ylab, col = "red")
    for(i in 1:length(alpha)) {
      pdf.x = dexp(x.rang, rate=1/theta[i])*alpha[i]
      lines(x.rang, pdf.x, lty = 2, col = "blue") }
  }
  
  if(is.matrix(x)) {
    xx = c()
    for(i in 1:dim(x)[1]) xx = c(xx, rep(x[i,1], x[i,2]))
    x = as.numeric(xx)
  }
  
  if(is.vector(x)) {
    zz = ss/max(x)
    hist(0.5+x*zz, nclass = nclass, probability = T,
         main=main, xlab = xlab, ylab = ylab)
    x.rang = seq(1, ss, 1)
    pdf.x = dmix.exp(x.rang, alpha, theta*zz)
    lines(x = pdf.x, col = "red")
    if(sub.pdf) {
      for(i in 1:length(alpha)) {
        pdf.x = dexp(x.rang, rate=1/theta[i]/zz)*alpha[i]
        lines(x = pdf.x, col = "blue")  } 
    } 
  }
}