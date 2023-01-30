plotmix.norm.new <-
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
      hist(x, freq= F, nclass = k, main=main, xlab = xlab, ylab = ylab)
      xx.grid = seq(min(x), max(x), (diff(range(yy))/20)/50)
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