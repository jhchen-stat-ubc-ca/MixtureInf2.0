#' plotmix.norm
#'
#' @description This function plots the mixture density, together with the histogram of the data
#' when the data generated from this density is given.
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
#' @examples n=200
#' alpha = c(.2, .5, .3)
#' mu = c(-1, 2, 5)
#' sigma = c(1, .5, 1.1)
#' x = rmix.norm(n, alpha, mu, sigma)
#' plotmix.norm(x, xx.grid = NULL, alpha, mu, sigma, m0=3,
#' k = 20, extra.height = 1.05, comp = TRUE, hist.ind = TRUE,
#' main="", xlab="Observed values", ylab="Density/Histogram")
#' @export
plotmix.norm <- function(x = NULL, xx.grid = NULL, alpha, mu, sigma, m0,
                         k = 20, extra.height = 1.05, comp = TRUE, hist.ind = TRUE,
                         main = "", xlab = "Observed values", ylab = "Density/Histogram") {
  if (hist.ind && is.null(x)) stop("data needed for histogram")
  if (hist.ind) {
    if (is.matrix(x)) x <- rep(x[, 1], x[, 2])
    if (is.vector(x)) {
      h <- hist(x, breaks = k, plot = FALSE)
      bin_width <- diff(h$breaks)[1]
      dens <- h$counts / (length(x) * bin_width)
      ylim_max <- max(dens) * extra.height
      hist(x, breaks = k, probability = TRUE, ylim = c(0, ylim_max),
           main = main, xlab = xlab, ylab = ylab)
      xx.grid <- seq(min(x), max(x), by = (diff(range(x)) / k) / 50)
      sig <- sqrt(sigma)
      sub.density <- t(sapply(1:m0, function(j) alpha[j] * dnorm(xx.grid, mu[j], sig[j])))
      mixture.density <- colSums(sub.density)
      lines(xx.grid, mixture.density, lty = 1)
      if (comp) for (j in 1:m0) lines(xx.grid, sub.density[j, ], lty = 2)
    }
  } else {
    sig <- sqrt(sigma)
    if (is.null(xx.grid)) {
      xx.min <- min(mu - 2.5 * sig)
      xx.max <- max(mu + 2.5 * sig)
      xx.grid <- seq(xx.min, xx.max, length.out = 1000)
    }
    sub.density <- t(sapply(1:m0, function(j) alpha[j] * dnorm(xx.grid, mu[j], sig[j])))
    mixture.density <- colSums(sub.density)
    plot(xx.grid, mixture.density, type = "l", main = main, xlab = xlab, ylab = ylab)
    if (comp) for (j in 1:m0) lines(xx.grid, sub.density[j, ], lty = 2)
  }
}
