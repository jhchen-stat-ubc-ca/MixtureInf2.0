### not meant for users, called by pmle.norm.R

pmle.norm.sub <-
function(x, m0, lambda, an, init.val, n.init, n.iter, max.iter, tol)
{
  # inputs are from pmle.norm, no need of checking.
		sn=var(x);  n=length(x)
		if (is.null(init.val)==F) n.init = 1
		
		output=c()
		for (i in 1:n.init) {
			if (is.null(init.val)) {
			    ## randomly generate initial values
				mu = sort(sample(x, m0))
				  tmp = (mu[-1] + mu[-m0])/2
				  alpha = rep(1, m0)
				  for(ii in 1:(m0-1)) alpha[ii] = sum(x<tmp[[ii]])/n
          alpha = alpha - c(0, alpha[-m0])
        sigma = c(mu, max(x))-c(min(x), mu)
				sigma = (sigma[-1]+sigma[-m0])/4
				sigma = sigma^2
			}	else {	### if initial mixing is provided
			  alpha = init.val[1:m0]
				mu = init.val[(m0+1):(2*m0)]
				sigma = init.val[(2*m0+1):(3*m0)]
			}
		  para0 = c(alpha, mu, sigma)
			for (j in 1:n.iter) {
			outpara = pmle.norm.sub.a(x, para0, lambda, an)
			para0 = outpara[1:3*m0]
			}  ###run n.iter EM-iterations first
			output=rbind(output,outpara[1:(3*m0)])	
		}
		index = which.max(output[,(3*m0+2)])
		para0 = output[index,1:(3*m0)]
		ploglike0 = output[index,(3*m0+2)]
		increment =1
		t=0
		### restart EM-iteration with the initial value 
		### that achieved the largest penalized log-likelihood
		while (increment > tol & t < max.iter) {
			outpara = pmle.norm.sub.a(x, para0, lambda, an)
			para0 = outpara[1:(3*m0)]
			ploglike1 = outpara[3*m0+2]
			increment = ploglike1 - ploglike0
			ploglike0 = ploglike1
			t = t+1
		}
		list('alpha'= outpara[1:m0],
		'means' = outpara[(m0+1):(2*m0)],
		'variances'= outpara[(2*m0+1):(3*m0)],
		'loglik'=outpara[3*m0+1],
		'ploglik'=outpara[3*m0+2])
}
