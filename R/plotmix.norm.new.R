plotmix.norm.new <-
  function(x = NULL, xx.grid = NULL, alpha, mu, sigma, m0,
           k = 20, extra.height = 1.05, comp = T, hist.ind = T,
           main="", xlab="Observed values", ylab="Density/Histogram") {
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
        sub.density = rbind(sub.density, alpha[j]*dnorm(xx.grid, mu[j], sig[j])) }
      mixture.density = colSums(sub.density)
      density.height = max(bb) * extra.height
      ### mixture and subpop densities.

      hist(x, prob=T, nclass = k, main=main, xlab = xlab, ylab = ylab)
      
      ### add densities
      lines(xx.grid, mixture.density, lty=1)
      if (comp) for (j in 1:m0) lines(xx.grid, sub.density[j,], lty=2)
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