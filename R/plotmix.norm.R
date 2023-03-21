#' plotmix.norm
#'
#' @param x The input data that can either be a vector or a matrix with the 1st column being the observed values
#'          and the 2nd column being the corresponding frequency.
#' @param xx.grid The grid of the histogram.
#' @param alpha A vector of the mixing proportions.
#' @param mu A vector of the component means.
#' @param sigma A vector of the component standard deviations.
#' @param m0 The order of the finite mixture model.
#' @param k The number of bars for the histogram.
#' @param extra.height Additional height multiplier used to enlarge the plot vertically.  
#' @param comp A parameter for the component fitted density. 
#'             comp=T means component fitted densities are drawn, and comp=F means no component fitted densities.
#' @param hist.ind It is true if there's data to be used to plot the histogram.
#' @param main The title of the graph.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#'
#' @return It returns the histogram of observations and the plot of the fitted density
#' @author Shaoting Li, Jiahua Chen and Pengfei Li
#'
#' @examples n=2000
#' alpha = c(.2, .5, .3)
#' mu = c(-1, 2, 5)
#' sigma = c(1, .5, 1.1)
#' yy = rmix.norm(n, alpha, mu, sigma)
#' plotmix.norm(yy, xx.grid = NULL, alpha, mu, sigma, m0=3,
#' k = 20, extra.height = 1.05, comp = T, hist.ind = T,
#' main="", xlab="Observed values", ylab="Density/Histogram")
plotmix.norm <-
  function(x = NULL, xx.grid = NULL, alpha, mu, sigma, m0,
           k = 20, extra.height = 1.05, comp = T, hist.ind = T,
           main="", xlab="Observed values", ylab="Density/Histogram") {
    if(hist.ind & is.null(x)) stop("data needed for histogram")
    
    if(hist.ind) {
      if(is.matrix(x)) {
        xx = c()
        for(i in 1:dim(x)[1]) xx = c(xx, rep(x[i,1], x[i,2]))
        x = as.numeric(xx)
      }
      if (is.vector(x)) {
        hist(x, freq= F, ylim = c(0,max((hist(x, freq = F, nclass = k))$density)*extra.height), nclass = k, main=main, xlab = xlab, ylab = ylab)
        xx.grid = seq(min(x), max(x), (diff(range(yy))/k)/50)
        sub.density = c(); sig = sigma^.5
        for (j in 1:m0) {
          sub.density = rbind(sub.density, alpha[j]*dnorm(xx.grid, mu[j], sig[j])) }
        mixture.density = colSums(sub.density) 
        lines(xx.grid, mixture.density, lty=1)
        if (comp) for (j in 1:m0) lines(xx.grid, sub.density[j,], lty=2)
      }
    }
    
    if(hist.ind==F) {
      sig = sigma^.5
      if(is.null(xx.grid)) {
        xx.min = min(mu - 2.5*sig[j]); xx.max = max(mu + 2.5*sig[j])
        xx.grid = seq(xx.min, xx.max, (xx.max-xx.min)/1000) }
      sub.density = c()
      for (j in 1:m0) {
        sub.density = rbind(sub.density, alpha[j]*dnorm(xx.grid, mu[j], sig[j])) }
      mixture.density = colSums(sub.density)
      ### mixture and subpop densities.
      plot(xx.grid, mixture.density, "l", main=main, xlab=xlab, ylab=ylab)
      for (j in 1:m0) lines(xx.grid, sub.density[j,], lty=2)
    }
  }