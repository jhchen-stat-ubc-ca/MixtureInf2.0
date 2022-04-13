penalty.norm.f <- function(alpha,sigma,sigma0,an,lambda)
{
	-an*sum(sigma0/sigma + log(sigma/sigma0))+lambda*sum(log(alpha))
}
