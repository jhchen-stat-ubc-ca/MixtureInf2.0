#' mom.bmm
#'
#' @description A sub-function for \code{pmle.beta}, performs the method of moments (MoM) estimation 
#'              for the beta mixture model.
#'
#'
#' @param x The input data, assumed to follow a beta mixture distribution.
#' @param m0 The number of components in the mixture.
#' @param maxit Maximum number of iterations allowed in the nonlinear solver.
#'
#' @return A matrix of estimated initial parameters with valid solutions.
#' @export
mom.bmm <- function(x, m0, maxit) {
  MM_BMM <- function(params) {
    m <- (length(params) + 1) / 3
    mix_porp <- params[1:(m - 1)]
    alpha <- params[m:(2 * m - 1)]
    beta  <- params[(2 * m):(3 * m - 1)]
    observed <- vapply(1:(3 * m - 1), function(k) mean(x^k), numeric(1))
    mm <- numeric(3 * m - 1)
    
    for (i in seq_along(mm)) {
      for (j in 1:(m - 1)) {
        mm[i] <- mm[i] + mix_porp[j] * prod((alpha[j] + 0:(i - 1)) / (alpha[j] + beta[j] + 0:(i - 1)))
      }
      mm[i] <- mm[i] + (1 - sum(mix_porp)) * prod((alpha[m] + 0:(i - 1)) / (alpha[m] + beta[m] + 0:(i - 1)))
    }
    (observed - mm)^2
  }
  
  repeat {
    kmeans_init <- ClusterR::KMeans_rcpp(matrix(x,ncol=1),m0)
    cluster_assignments <- kmeans_init$clusters
    mix_porp <- table(cluster_assignments) / length(cluster_assignments)
    
    theta <- mom.calculation(x, cluster_assignments, m0)
    
    test_out <- nleqslv::testnslv(c(mix_porp[1:(m0 - 1)], theta), MM_BMM, control = list(maxit = maxit))$out
    valid <- test_out$termcd < 4
    if (!any(valid)) next
    
    methods <- test_out$Method[valid]
    globals <- test_out$Global[valid]
    
    solns <- lapply(seq_along(methods), function(i) {
      nleqslv::nleqslv(c(mix_porp[1:(m0 - 1)], theta), MM_BMM,
                       method = methods[i], global = globals[i],
                       control = list(maxit = maxit))
    })
    
    valid_solns <- vapply(solns, function(sol) {
      sum(sol$x[1:(m0 - 1)]) < 1 && all(sol$x > 0)
    }, logical(1))
    
    if (any(valid_solns)) {
      return(do.call(cbind, solns[valid_solns]))
    }
  }
}


#' mom.calculation
#'
#' @description A sub-function for \code{mom.bmm}, generates starting values for the method 
#'              of moments estimation of a beta mixture. 
#'
#' @param x The input data, assumed to be from a beta mixture.
#' @param cluster_assignments A vector of component identity that identifies which subpopulation that this observed value belongs to.
#' @param m0 The number of components in the mixture.
#'
#' @return A numeric vector containing concatenated alpha and beta estimates for each component.
#' @export

mom.calculation <- function(x, cluster_assignments, m0) {
  alpha <- numeric(m0)
  beta <- numeric(m0)
  
  for (i in 1:m0) {
    cluster_data <- x[cluster_assignments == i]
    temp <- mom.beta(cluster_data)
    alpha[i] <- temp$alpha
    beta[i]  <- temp$beta
  }
  
  c(alpha, beta)
}

#' mom.beta
#'
#' @description A sub-function for \code{mom.calculation}, estimates the parameters of a single beta 
#'              distribution using the method of moments.
#'
#'
#' @param data A numeric vector assumed to be sampled from a beta distribution.
#'
#' @return A list with components:
#' \item{alpha}{The estimated alpha parameter.}
#' \item{beta}{The estimated beta parameter.}
#' @export

mom.beta <- function(data) {
  mu <- mean(data)
  var_data <- var(data)
  alpha <- mu * (mu * (1 - mu) / var_data - 1)
  beta  <- (1 - mu) * (mu * (1 - mu) / var_data - 1)
  return(list(alpha = alpha, beta = beta))
}
