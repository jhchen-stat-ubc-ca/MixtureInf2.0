#' emstat.exp
#'
#' @description This function computes the EM-test statistics for the exponential mixture.
#' It is used in the emtest.exp function.
#' 
#' @export
emstat.exp <- function(x, outnull, CC, n.init, n.iter, tol, k, 
                       max.iter) {
  theta0 = outnull$theta
  m0     = length(theta0)	
  ### create the collection of beta	
  bbeta=c()
  for(h in 1:m0)  {
    bbeta = rbind(cbind(bbeta, rep(0.1,3^{h-1})),
                  cbind(bbeta, rep(0.3,3^{h-1})),
                  cbind(bbeta, rep(0.5,3^{h-1})))
  }
  
  pln1 = para = c()
  ##For each beta, calculate the statistic m_n^{(k)}
  for(j in 1:(3^m0)) {
    output=c()
    beta.j = bbeta[j,]   ### avoid R-function beta.
    para0 = maxmm.exp(x, beta.j, theta0, n.init, n.iter, tol, max.iter)
    alpha  = para0$alpha
    theta1 = para0$theta1
    theta2 = para0$theta2
    
    ### EM-iteration for EM-test
    for(i in 1:k) {
      alpha1 = alpha*beta.j	
      alpha2 = alpha*(1-beta.j)
      pdf.part1 = apply(as.matrix(1/theta1,ncol=1),1,dexp,x=x)
      pdf.part2 = apply(as.matrix(1/theta2,ncol=1),1,dexp,x=x)
      pdf.sub1 = t(t(pdf.part1)*alpha1)+1e-100/m0
      pdf.sub2 = t(t(pdf.part2)*alpha2)+1e-100/m0
      pdf = apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)	
      w1 = pdf.sub1/pdf	
      w2 = pdf.sub2/pdf
      
      ###M-step of EM-algorithm
      for(h in 1:m0)  {
        if (sum(w1[,h])/(sum(w1[,h])+sum(w2[,h]))<=0.5)
          beta.j[h]=min((sum(w1[,h])+CC)/(sum(w1[,h])+sum(w2[,h])+CC),0.5)
        else
          beta.j[h]=max((sum(w1[,h]))/(sum(w1[,h])+sum(w2[,h])+CC),0.5)
      }
      alpha=apply(w1+w2,2,mean)
      alpha1=alpha*beta.j
      alpha2=alpha*(1-beta.j)
      theta1=apply(w1*x,2,sum)/apply(w1,2,sum)
      theta2=apply(w2*x,2,sum)/apply(w2,2,sum)
    }		
    ###compute the penalized log-likelihood value and EM-test statistic
    pdf.part1=apply(as.matrix(1/theta1,ncol=1),1,dexp,x=x)
    pdf.part2=apply(as.matrix(1/theta2,ncol=1),1,dexp,x=x)
    pdf.sub1=t(t(pdf.part1)*alpha1)+1e-100/m0	
    pdf.sub2=t(t(pdf.part2)*alpha2)+1e-100/m0
    pdf=apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)
    para=rbind(para,c(alpha1,alpha2,theta1,theta2))
    pln1[j] = sum(log(pdf)) + CC*sum(log(1-abs(1-2*beta.j)))
  }
  
  index = which.max(pln1)  ## pick the winner.
  para1 = para[index,]
  emnk  = 2*(pln1[index]-outnull$loglik)
  c(emnk, para1)           ## report the EM-stat + parameter.
}