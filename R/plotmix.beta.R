#' plotmix.beta
#'
#' @description This function plots the mixture density, together with the histogram of the data
#' when the data generated from this density is given.
#' @param x The input data which is a vector.
#' @param xx.grid The grid of the histogram.
#' @param porp A vector of the mixing proportions.
#' @param alpha A vector of the component alphas.
#' @param beta A vector of the component betas.
#' @param m0 The order of the finite mixture model.
#' @param k The number of bars for the histogram.
#' @param extra.height Additional height multiplier used to enlarge the plot vertically.  
#' @param comp A parameter for the component fitted density. 
#'             comp=T means component fitted densities are drawn, and comp=F means no component fitted densities.
#' @param main The title of the graph.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#'
#' @return It returns the histogram of observations and the plot of the fitted density
#'
#' @examples n=2000
#' porp = c(.2, .5, .3)
#' alpha = c(3, 2, 5)
#' beta = c(1, .5, 1.1)
#' x = rmix.beta(n, porp, alpha, beta)
#' plotmix.beta(x, xx.grid=NULL, porp, alpha, beta, m0=3,
#' main="", xlab="Observed values", ylab="Density/Histogram")
#' @export
plotmix.beta <-
  function(x, xx.grid = NULL, porp, alpha, beta, m0,
           k = 20, extra.height = 1.05, comp = T,
           main="", xlab="Observed values", ylab="Density/Histogram") {
    if(is.matrix(x)) {
      xx = c()
      for(i in 1:dim(x)[1]) xx = c(xx, rep(x[i,1], x[i,2]))
      x = as.numeric(xx)
    }
    if (is.vector(x)) {
      hist(x, freq= F, ylim = c(0,max((hist(x, freq = F, nclass = k))$density)*extra.height), 
           nclass = k, main=main, xlab = xlab, ylab = ylab)
      xx.grid = seq(min(x), max(x), (diff(range(x))/k)/50)
      sub.density = c()
      for (j in 1:m0) {
        sub.density = rbind(sub.density, porp[j]*dbeta(xx.grid, alpha[j], beta[j])) }
      mixture.density = colSums(sub.density)
      lines(xx.grid, mixture.density,lwd=2)
      if (comp==T) { for (j in 1:m0) 
        lines(xx.grid, sub.density[j,],lwd=1,col=j+1)}
    }
  }