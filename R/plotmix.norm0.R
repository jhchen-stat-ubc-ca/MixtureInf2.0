#' plotmix.norm0
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
#' @examples n=3000
#' alpha = c(.2, .4, .4)
#' mu = c(-1, 2, 5)
#' sigma = c(1, .5, 1.1)
#' yy = rmix.norm(n, alpha, mu, sigma)
#' plotmix.norm0(yy, xx.grid = NULL, alpha, mu, sigma, m0=3,
#' k = 20, extra.height = 1.05, comp = T, hist.ind = T,
#' main="", xlab="Observed values", ylab="Density/Histogram")
plotmix.norm0 <- function(x = NULL, xx.grid = NULL, alpha, mu, sigma, m0,
                          k = 20, extra.height = 1.05, comp = T, hist.ind = T,
                          main="", xlab="Observed values", ylab="Density/Histogram")
{ 
  if(hist.ind & is.null(x)) stop("data needed for histogram")
  
  if(hist.ind) {
    if (is.vector(x)) { 
      xx = x
      nn = length(xx);   xx = sort(xx)
      rr = (xx[nn]-xx[1])/k
      aa = xx[1] + (0:k)*rr
      bb = rep(0, k+1)
      for (i in 1:(k+1)) bb[i] = sum(xx >= aa[i]-rr/2 & xx < aa[i]+rr/2)
      bb = bb/nn/rr
      xx.grid = seq(aa[1]-rr , aa[k+1]+rr, rr/50) 
    }  ## aa: grid for histogram; bb: height, rr: width
    
    if (is.matrix(x)) {  ## handle frequency data
      k = dim(x)[1]-1
      aa = x[,1] 
      rr = (aa[2]-aa[1])
      bb = x[,2]/sum(x[,2])/rr
      xx.grid = seq(aa[1]-rr , aa[k+1]+rr, rr/50)
    }
    
    sub.density = c(); sig = sigma^.5
    for (j in 1:m0) {
      sub.density = rbind(sub.density, alpha[j]*dnorm(xx.grid, mu[j], sig)) }
    mixture.density = colSums(sub.density)
    density.height = max(bb) * extra.height
    ### mixture and subpop densities.
    
    plot(c(aa[1]-rr, aa[k+1]+rr), c(0, density.height), "n", 
         main=main, xlab=xlab, ylab=ylab)
    lines(c(aa[1]-rr,aa[k+1]+rr), c(0,0), "l", col="black")
    lines(aa[1]+c(-0.5, -0.5)*rr, c(0,  bb[1]),"l",col="blue")
    for (i in 1:k) {
      lines(aa[i]+c(-0.5, 0.5)*rr, c(bb[i],  bb[i]),"l",col="blue")
      lines(aa[i]+c(0.5, 0.5)*rr, c(bb[i], bb[i+1]),"l",col="blue") }
    lines(aa[k+1]+c(-0.5, 0.5)*rr, c(bb[k+1],  bb[k+1]),"l",col="blue")
    lines(aa[k+1]+c(0.5, 0.5)*rr, c(bb[k+1],0),"l",col="blue")
    
    ### add densities
    lines(xx.grid, mixture.density, lty=1)
    if (comp) for (j in 1:m0) lines(xx.grid, sub.density[j,], lty=2)
  }
  
  if(hist.ind==F) {
    sig = sigma^.5
    if(is.null(xx.grid)) {
      xx.min = min(mu - 2.5*sig); xx.max = max(mu + 2.5*sig)
      xx.grid = seq(xx.min, xx.max, (xx.max-xx.min)/1000) }
    sub.density = c()
    for (j in 1:m0) {
      sub.density = rbind(sub.density, alpha[j]*dnorm(xx.grid, mu[j], sig)) }
    mixture.density = colSums(sub.density)
    ### mixture and subpop densities.
    plot(xx.grid, mixture.density, "l", main=main, xlab=xlab, ylab=ylab)
    for (j in 1:m0) lines(xx.grid, sub.density[j,], lty=2)
  }
}