#' emstat.norm0
#'
#' @description The function computes the EM-test statistics for univariate normal mixture
#' with known variance. It is used in the emtest.norm0 function.
emstat.norm0 <-
  function(xx, mle0, ln0, m0, CC, n.init, n.iter, tol, k=3)
    # xx:        data as a plain vector. 
    # mle0:     output from pmle.norm0 function without penalty. 
    ##  Two rows: mixing proportions and subpop means.
    # ln0:    penalized log likelihood at fitted null.
    # m0:      the null order
    # CC:        tuning parameter for EM-test procedure. 
    # n.init, n.iter: for EMtest.maxmm.norm0.
    # tol:     tolerance for increment of log likelihood
    # k:	       number of iterations for the EM-test statistic.
  {
    bbeta = c()
    for(ii in 1:m0) {
      bbeta = rbind(cbind(bbeta, rep(0.1,3^{ii-1})),
                    cbind(bbeta, rep(0.3,3^{ii-1})),
                    cbind(bbeta, rep(0.5,3^{ii-1})))  }
    ### all possible mixing proportions combinations from
    ###    {0.1, 0.3, 0.5}  Li and Chen (JASA 2010)
    ###   Total is 3^m0.
    pln1 = para = c()
    ##For each beta, calculate the statistic M_n^{(k)}
    for(j in 1:(3^m0)) {
      para0 = EMtest.maxmm.norm0(xx, bbeta[j,], mle0, m0, 
                                 n.init, n.iter, tol)
      ### pmle of mixture with special strcture. (JASA1085 R)
      theta1=para0$theta1
      theta2=para0$theta2
      alpha1= para0$alpha*bbeta[j,]	
      alpha2= para0$alpha*(1 - bbeta[j,])
      
      ### K EM-iterations as required by EM-test
      for(i in 1:k) {
        ###E-step
        pdf.part1=apply(as.matrix(theta1,ncol=1),1,dnorm,x=xx)
        pdf.part2=apply(as.matrix(theta2,ncol=1),1,dnorm,x=xx)
        pdf.sub1 = t(t(pdf.part1)*alpha1)+1e-100/m0
        pdf.sub2 = t(t(pdf.part2)*alpha2)+1e-100/m0
        pdf = apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)	
        w1  = pdf.sub1/pdf; 	w2 = pdf.sub2/pdf
        
        ###M-step
        for(h in 1:m0) {
          tt1 = sum(w1[,h]); tt2 = sum(w2[,h])
          if (tt1 <= tt2) {
            bbeta[j,h]=min((tt1 + CC)/(tt1+tt2+CC),0.5)
          } else {
            bbeta[j,h]=max(tt1/(tt1+tt2+CC),0.5) }
        }   ## loop h
        alpha = apply(w1+w2, 2, mean)
        alpha1 = alpha*bbeta[j,];  alpha2 = alpha*(1 - bbeta[j,])
        theta1 = apply(w1*xx,2,sum)/apply(w1,2,sum)
        theta2 = apply(w2*xx,2,sum)/apply(w2,2,sum)
      }  ## loop i
      
      ###Compute the penalized log-likelihood and other values
      pdf.part1=apply(as.matrix(theta1,ncol=1),1,dnorm,x=xx)
      pdf.part2=apply(as.matrix(theta2,ncol=1),1,dnorm,x=xx)
      pdf.sub1 =t(t(pdf.part1)*alpha1)+1e-100/m0	
      pdf.sub2 =t(t(pdf.part2)*alpha2)+1e-100/m0
      pdf = apply(pdf.sub1,1,sum)+apply(pdf.sub2,1,sum)
      para = rbind(para,c(alpha1,alpha2,theta1,theta2))
      pln1 = c(pln1, sum(log(pdf))+CC*sum(log(1-abs(1-2*bbeta[j,]))))
    } ## loop j.
    
    index = which.max(pln1)
    para1 = para[index,]
    emnk = 2*(pln1[index] - ln0)
    
    c(emnk, para1)
  }