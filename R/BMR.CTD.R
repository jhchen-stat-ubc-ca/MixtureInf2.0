#' beta.Hellinger
#'
#' @description Computes the squared Hellinger distance between two beta distributions with 
#' parameters \code{(a1, b1)} and \code{(a2, b2)}. It is a sub function used in the beta.barycenter function.
#'
#'
#' @param a1 First shape parameter of the first beta distribution.
#' @param b1 Second shape parameter of the first beta distribution.
#' @param a2 First shape parameter of the second beta distribution.
#' @param b2 Second shape parameter of the second beta distribution.
#'
#' @return A numeric value representing the squared Hellinger distance.
#' @export
beta.Hellinger <- function(a1, b1, a2, b2) {
  return(1 - beta((a1 + a2) / 2, (b1 + b2) / 2) / sqrt(beta(a1, b1) * beta(a2, b2)))  
}

#' beta.SKL
#'
#' @description Computes the symmetrized KL divergence between two beta distributions with 
#' parameters \code{(a1, b1)} and \code{(a2, b2)}. It is a sub function used in the beta.barycenter function.
#'
#'
#' @param a1 First shape parameter of the first beta distribution.
#' @param b1 Second shape parameter of the first beta distribution.
#' @param a2 First shape parameter of the second beta distribution.
#' @param b2 Second shape parameter of the second beta distribution.
#'
#' @return A numeric value representing the symmetrized KL divergence.
#' @export
beta.SKL <- function(a1, b1, a2, b2) {
  kl12 <- lbeta(a2,b2) - lbeta(a1,b1) +
    (a1-a2)*digamma(a1) + (b1-b2)*digamma(b1) -
    ((a1-a2)+(b1-b2))*digamma(a1+b1)
  kl21 <- lbeta(a1,b1) - lbeta(a2,b2) +
    (a2-a1)*digamma(a2) + (b2-b1)*digamma(b2) -
    ((a2-a1)+(b2-b1))*digamma(a2+b2)
  return(kl12 + kl21)
}

#' beta.barycenter
#'
#' @description Computes the barycenter of a collection of beta distributions 
#' with respect to a divergence—the Hellinger distance or the symmetric KL divergence. This function is used in \code{BMR.CTD}.
#'
#' @param assignments An integer vector indicating the cluster assignment of each original beta component.
#' @param M The number of target components in the reduced mixture.
#' @param cluster_weights A numeric vector of weights for each original beta component.
#' @param cluster_alphas A numeric vector of \eqn{\alpha} parameters for each original beta component.
#' @param cluster_betas A numeric vector of \eqn{\beta} parameters for each original beta component.
#' @param divergence Character string, either "Hellinger" (default) or "SKL"
#'
#' @return A matrix with \code{M} rows and 2 columns, where each row contains the estimated 
#'         \code{(alpha, beta)} parameters of the corresponding cluster's barycenter.
#' @export
beta.barycenter <- function(assignments, M, cluster_weights,
                            cluster_alphas, cluster_betas,
                            divergence = "Hellinger") {
  dist_fun <- if (divergence == "Hellinger") beta.Hellinger else beta.SKL
  res <- matrix(NA_real_, nrow = M, ncol = 2)
  for (m in seq_len(M)) {
    idx <- which(assignments == m)
    if (!length(idx)) next
    w <- cluster_weights[idx]
    a <- cluster_alphas[idx]
    b <- cluster_betas[idx]
    mu <- sum(w * a / (a + b))
    start <- c(max(mu, 1e-3), max(1 - mu, 1e-3))
    opt <- optim(
      par = start,
      fn = function(par) sum(w * dist_fun(a, b, par[1], par[2])),
      method = "L-BFGS-B",
      lower = c(1e-6, 1e-6),
      control = list(maxit = 1000, factr = 1e6)
    )
    res[m, ] <- opt$par
  }
  res
}


#' init_reduced_beta_mixture
#'
#' @description Generates initial parameters for a reduced beta mixture by sampling from the original 
#' mixture and clustering the samples with k-means. The Score-Adjusted Moment estimation of the clustered samples 
#' are then used to estimate the parameters. It is a sub function used in the BMR.CTD function.
#'
#' @param weights A numeric vector of mixture weights.
#' @param alphas A numeric vector of \eqn{\alpha} parameters for each Beta component.
#' @param betas A numeric vector of \eqn{\beta} parameters for each Beta component.
#' @param M The number of components in the reduced mixture.
#' @param n_sample The number of samples to draw from the original mixture. Default is 10000.
#'
#' @return A list containing the initialized \code{alphas}, \code{betas}, and \code{weights} 
#'         for the reduced Beta mixture.
#' @export
init_reduced_beta_mixture <- function(weights, alphas, betas, M, n_sample = 10000) {
  samples <- rmix.beta(n_sample, weights, alphas, betas)
  km <- ClusterR::KMeans_rcpp(matrix(samples, ncol = 1), clusters = M, initializer = "kmeans++")
  ids <- km$clusters
  props <- tabulate(ids, nbins = M) / n_sample
  pars <- vapply(seq_len(M), function(m) {
    xj <- samples[ids == m]
    est <- sam.beta(xj)
    c(est$alpha, est$beta)
  }, numeric(2))
  
  list(weights = props, alphas = pars[1, ], betas = pars[2, ])
}




#' BMR.CTD
#'
#' @description Performs beta mixture reduction (BMR) via Composite Transportation Divergence (CTD),
#' which iteratively reassigns components based on a divergence (the Hellinger distance or symmetric KL divergence), and recomputes cluster barycenters 
#' to minimize total cost.
#'
#' @param orig_weights A numeric vector of original mixture weights.
#' @param orig_alphas A numeric vector of \eqn{\alpha} parameters for the original beta components.
#' @param orig_betas A numeric vector of \eqn{\beta} parameters for the original beta components.
#' @param M The number of components in the reduced mixture.
#' @param divergence Character string, either "Hellinger" (default) or "SKL"
#' @param max_iter Maximum number of iterations for the reduction algorithm. Default is 100.
#' @param tol Tolerance threshold for convergence based on parameter change. Default is \code{1e-6}.
#' @param n_sample Number of samples used for initial reduction parameter estimation. Default is 10000.
#'
#' @return A list with the following components:
#' \item{reduced_weights}{Estimated weights of the reduced components.}
#' \item{reduced_alphas}{Estimated \eqn{\alpha} parameters for the reduced mixture.}
#' \item{reduced_betas}{Estimated \eqn{\beta} parameters for the reduced mixture.}
#' \item{total_cost}{Total cost of the final clustering solution.}
#' \item{assignments}{Cluster assignments of original components.}
#' \item{n_iter}{Number of iterations until convergence.}
#' \item{converged}{Logical indicating whether convergence was achieved.}
#' \item{divergence_used}{The choice of the distance/divergence used}
#'
#' @examples
#' orig_weights <- c(0.3, 0.3, 0.3, 0.1)
#' orig_alphas <- c(3, 4, 3.5, 7)
#' orig_betas <- c(5, 6, 5.5, 2)
#' result <- BMR.CTD(orig_weights, orig_alphas, orig_betas, M = 2)
#' curve({dmix.beta(x, orig_weights, orig_alphas, orig_betas)}, from = 0, to = 1,
#'       col = "blue", lwd = 2, xlab = "x", ylab = "Density", main = "Original vs Reduced Mixture")
#' curve({dmix.beta(x, result$mix_prop, result$alpha, result$beta)},
#'       add = TRUE, col = "red", lwd = 2, lty = 2)
#' @export
BMR.CTD <- function(orig_weights, orig_alphas, orig_betas, M,
                    divergence = "Hellinger",
                    max_iter = 100, tol = 1e-6, n_sample = 10000) {
  
  dist_fun <- if (divergence == "Hellinger") beta.Hellinger else beta.SKL
  orig_weights <- orig_weights / sum(orig_weights)
  
  make_cost <- function(ra, rb) {
    m <- length(ra)
    n <- length(orig_alphas)
    cost <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) {
      cost[, j] <- dist_fun(orig_alphas, orig_betas, ra[j], rb[j])
    }
    cost
  }
  
  init <- init_reduced_beta_mixture(orig_weights, orig_alphas, orig_betas, M, n_sample)
  p_old <- c(init$weights, init$alphas, init$betas)
  
  for (iter in seq_len(max_iter)) {
    mw <- p_old[1:M]
    ra <- p_old[(M + 1):(2 * M)]
    rb <- p_old[(2 * M + 1):(3 * M)]
    
    cost <- make_cost(ra, rb)
    assign <- max.col(-cost)
    
    new_rw <- vapply(seq_len(M), function(m) sum(orig_weights[assign == m]), numeric(1))
    if (sum(new_rw) == 0) new_rw[] <- 1 / M
    new_rw <- new_rw / sum(new_rw)
    
    new_ab <- beta.barycenter(assign, M, orig_weights, orig_alphas, orig_betas, divergence)
    new_ra <- new_ab[, 1]
    new_rb <- new_ab[, 2]
    
    p_new <- c(new_rw, new_ra, new_rb)
    if (anyNA(p_new) || max(abs(p_new - p_old)) < tol) break
    p_old <- p_new
  }
  
  mw <- p_old[1:M]
  ra <- p_old[(M + 1):(2 * M)]
  rb <- p_old[(2 * M + 1):(3 * M)]
  mw <- mw / sum(mw)
  
  cost <- outer(seq_along(orig_alphas), seq_len(M),
                Vectorize(function(i, j)
                  dist_fun(orig_alphas[i], orig_betas[i], ra[j], rb[j])))
  assign <- max.col(-cost)
  total_cost <- sum(orig_weights * cost[cbind(seq_along(assign), assign)])
  
  list(
    mix_prop = rousignif(mw),
    alpha = rousignif(ra),
    beta = rousignif(rb),
    total_cost = total_cost,
    assignments = assign,
    n_iter = iter,
    converged = (iter < max_iter),
    divergence_used = divergence
  )
}
