#' EMstat.pois
#'
#' @description This function computes the EM-test statistics for Poisson mixture. 
#' It is used in the emtest.pois function.
EMstat.pois <- function(count, freq, alpha0, theta0, m0, ln0, 
                        CC, n.init, n.iter, tol, k, max.iter)
{
  ### all possible mixing proportions combinations from
  ###    {0.1, 0.3, 0.5}  Li and Chen (JASA 2010)
  ###   Total is 3^m0.
  bbeta=c()
  for(h in 1:m0)  {
    bbeta=rbind(cbind(bbeta,rep(0.1,3^{h-1})),
                cbind(bbeta, rep(0.3,3^{h-1})),
                cbind(bbeta, rep(0.5,3^{h-1})))  }
  
  pln1 = para=c()
  ##For each beta, calculate the statistic M_n^{(k)}
  for(j in 1:(3^m0)) {
    output=c()
    beta.i = bbeta[j,]
    para0 = EMtest.maxmm.pois(count, freq, beta.i, 
                              alpha0, theta0, m0, n.init, n.iter, tol, max.iter)
    alpha1  = para0$alpha1
    alpha2  = para0$alpha2
    theta1 = para0$theta1
    theta2 = para0$theta2
    
    ### Iteration k times
    for(i in 1:k)	{
      pdf.part1=apply(as.matrix(theta1,ncol=1),1,dpois,x=count)
      pdf.part2=apply(as.matrix(theta2,ncol=1),1,dpois,x=count)
      pdf.component1=t(t(pdf.part1)*alpha1)+1e-100/m0
      pdf.component2=t(t(pdf.part2)*alpha2)+1e-100/m0
      pdf=apply(pdf.component1,1,sum)+apply(pdf.component2,1,sum)	
      w1=pdf.component1/pdf	
      w2=pdf.component2/pdf
      for(h in 1:m0) {
        tt1 = sum(freq*w1[,h]);  tt2 = sum(freq*w2[,h])
        if (tt1 < tt2) { beta.i[h]=min( (tt1+CC)/(tt1+tt2+CC), 0.5)
        } else { beta.i[h]=max(tt1/(tt1+tt2+CC), 0.5) }
      }
      alpha = apply(freq*w1+freq*w2,2,sum)/sum(freq)
      alpha1= alpha*beta.i
      alpha2= alpha*(1 - beta.i)
      theta1= apply(freq*w1*count,2,sum)/apply(freq*w1,2,sum)
      theta2= apply(freq*w2*count,2,sum)/apply(freq*w2,2,sum)
    }	 ### complete k iterations here.
    
    ###Compute the penalized log-likelihood value
    pdf = dmix.pois(count, c(alpha1, alpha2), c(theta1, theta2))
    para = rbind(para, c(alpha1,alpha2,theta1,theta2))
    pln1[j] = sum(freq*log(pdf)) + CC*sum(log(1-abs(1-2*beta.i)))
  }
  
  ## pick the winner
  index = which.max(pln1)
  para1 = para[index,]
  emnk = 2*(pln1[index] - ln0)
  
  c(emnk, para1)
}