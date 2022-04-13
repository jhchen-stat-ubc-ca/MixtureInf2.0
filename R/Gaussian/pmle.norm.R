### This function computes penalized MLE
###     under uni-variate finite Gaussian Mixture.
###   It organize input for pmle.norm.sub.R
pmle.norm <-
function(x, m0 = 2, lambda = 1, an = NULL, init.val = NULL,
         n.init=NULL, n.iter=50, tol=1e-8, rformat = FALSE)
{
  # x: 		data, either a vector or a matrix with two columns.
  #            col 1: observed, col 2: frequency 
  # m0:	 the order of mixture to be fitted.
  # lambda: modification of chen, chen, Kalbfleish
  # an: penalty on variance of Chen, Tan and Zhang; recommended n^{-1/2}
  # init.val:	user provided initial mixing distribution as a 3 X m0 matrix
  #          including mixing probability, component means and variances.
  # n.init:	user may us randomly generate n.init initials for the EM-algorithm.	
  # n.iter: the least number of EM iterations for each initial values.
  # tol:		EM concludes if the increment in p-likelihood is below.
  # rformat:	format for output, 
      # rformat=T: default of R-cran.
      #	rformat=F: if output > 0.001, report round(output,3);
      #            otherwise, report signif(output,3).
  if(m0==1) stop("this function is for mixture of order 2 or higher")
  
	if (is.data.frame(x)) {	
		if (ncol(x)==2) x=as.matrix(x)
		if (ncol(x)==1 | ncol(x)>2) x=x[,1]
	}
  
	if (is.matrix(x)) {
		xx=c()
		for (i in 1:nrow(x)) xx=c(xx, rep(x[i,1],x[i,2]))
	}       ## not ideal when data are discretized. 
          ## it expands it into a vector.
  if(is.null(an)) an = length(xx)^(0.5)
          ## default penalty size n^{-1/2}
	out = pmle.norm.sub(xx, m0, lambda, an, init.val, n.init, n.iter, tol)
	    ## leave the computation to the other function.
	alpha = out[[1]]
	mean.mix = out[[2]]
	var.mix = out[[3]]
	loglik = out[[4]]
	ploglik = out[[5]]
	
	if (rformat==F) {
		alpha.mix = rousignif(alpha)
		mean.mix = rousignif(mean.mix)
		var.mix = rousignif(var.mix)
		loglik = rousignif(loglik)
		ploglik = rousignif(ploglik)
	}

		list('PMLE of mixing proportions' = alpha,
		  'PMLE of means' = mean.mix,
		  'PMLE of variances' = var.mix,
		  'log-likelihood' = loglik,
		  'Penalized log-likelihood'= ploglik)
}
### Looks fine but need final approval.