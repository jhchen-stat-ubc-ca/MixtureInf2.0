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

#' beta.JS
#'
#' @description Computes the Jensen–Shannon divergence between two beta distributions with 
#' parameters \code{(a1, b1)} and \code{(a2, b2)}. It is a sub function used in the beta.barycenter function.
#'
#'
#' @param a1 First shape parameter of the first beta distribution.
#' @param b1 Second shape parameter of the first beta distribution.
#' @param a2 First shape parameter of the second beta distribution.
#' @param b2 Second shape parameter of the second beta distribution.
#'
#' @return A numeric value representing the Jensen–Shannon divergence.
#' @export
beta.JS <- function(a1, b1, a2, b2) {
  kl12 <- lbeta(a2,b2) - lbeta(a1,b1) +
    (a1-a2)*digamma(a1) + (b1-b2)*digamma(b1) -
    ((a1-a2)+(b1-b2))*digamma(a1+b1)
  kl21 <- lbeta(a1,b1) - lbeta(a2,b2) +
    (a2-a1)*digamma(a2) + (b2-b1)*digamma(b2) -
    ((a2-a1)+(b2-b1))*digamma(a2+b2)
  return(0.5 * (kl12 + kl21))
}

#' beta.barycenter
#'
#' @description Computes the barycenter of a collection of beta distributions 
#' with respect to a hybrid divergence—a convex combination of the squared Hellinger distance 
#' and the Jensen–Shannon divergence—using a gradient-free optimization method.
#' This function is a subroutine used in \code{BMR.CTD}.
#'
#' @param assignments An integer vector indicating the cluster assignment of each original beta component.
#' @param M The number of target clusters (components) in the reduced mixture.
#' @param cluster_weights A numeric vector of weights for each original beta component.
#' @param cluster_alphas A numeric vector of \eqn{\alpha} parameters for each original beta component.
#' @param cluster_betas A numeric vector of \eqn{\beta} parameters for each original beta component.
#' @param zeta A scalar between 0 and 1 representing the weight on the squared Hellinger distance 
#'        (with \code{1 - zeta} on the Jensen–Shannon divergence).
#'
#' @return A matrix with \code{M} rows and 2 columns, where each row contains the estimated 
#'         \code{(alpha, beta)} parameters of the corresponding cluster's barycenter.
#' @export
beta.barycenter <- function(assignments, M, cluster_weights, cluster_alphas, 
                            cluster_betas, zeta) {
  result <- matrix(NA_real_, nrow = M, ncol = 2)
  
  for (m in seq_len(M)) {
    idx <- which(assignments == m)
    if (length(idx) == 0) next
    
    w <- cluster_weights[idx]
    a <- cluster_alphas[idx]
    b <- cluster_betas[idx]
    
    mu0 <- sum(w * (a / (a + b)))
    init_par <- c(mu0, 1 - mu0)
    
    res <- optim(
      par = init_par,
      fn = function(par) {
        hd  <- beta.Hellinger(a, b, par[1], par[2])
        jsd <- beta.JS       (a, b, par[1], par[2])
        sum(w * (zeta * hd + (1 - zeta) * jsd))
      },
      method = "L-BFGS-B",
      lower = c(1e-6, 1e-6),
      control = list(maxit = 1000, factr = 1e6)
    )
    
    result[m, ] <- res$par
  }
  result
}


#' init_reduced_beta_mixture
#'
#' @description Generates initial parameters for a reduced beta mixture by sampling from the original 
#' mixture and clustering the samples with k-means. The moments of the clustered samples 
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
  km <- kmeans(samples, centers = M)
  
  cluster_ids <- km$cluster
  proportions <- tabulate(cluster_ids, nbins = M) / n_sample
  
  init_params <- vapply(seq_len(M), function(m) {
    cluster_samples <- samples[cluster_ids == m]
    mu <- mean(cluster_samples)
    v  <- var(cluster_samples)
    s <- mu * (1 - mu) / v - 1
    alpha <- max(mu * s, 1e-3)
    beta  <- max((1 - mu) * s, 1e-3)
    
    c(alpha, beta)
  }, numeric(2))
  
  list(
    alphas  = init_params[1, ],
    betas   = init_params[2, ],
    weights = proportions
  )
}


#' BMR.CTD
#'
#' @description Performs beta mixture reduction (BMR) via Composite Transportation Divergence (CTD),
#' which iteratively reassigns components based on a hybrid divergence (a convex combination of the 
#' squared Hellinger distance and Jensen–Shannon divergence), and recomputes cluster barycenters 
#' to minimize total cost.
#'
#' @param orig_weights A numeric vector of original mixture weights.
#' @param orig_alphas A numeric vector of \eqn{\alpha} parameters for the original beta components.
#' @param orig_betas A numeric vector of \eqn{\beta} parameters for the original beta components.
#' @param M The number of components in the reduced mixture.
#' @param zeta Optional scalar between 0 and 1 specifying the relative weight on the Hellinger distance.
#'        If \code{NULL}, the optimal \code{zeta} is automatically selected to minimize the total cost.
#' @param max_iter Maximum number of iterations for the reduction algorithm. Default is 100.
#' @param tol Tolerance threshold for convergence based on parameter change. Default is \code{1e-6}.
#' @param n_sample Number of samples used for initial reduction parameter estimation. Default is 10000.
#'
#' @return A list with the following components:
#' \item{reduced_alphas}{Estimated \eqn{\alpha} parameters for the reduced mixture.}
#' \item{reduced_betas}{Estimated \eqn{\beta} parameters for the reduced mixture.}
#' \item{reduced_weights}{Estimated weights of the reduced components.}
#' \item{total_cost}{Total cost of the final clustering solution.}
#' \item{assignments}{Cluster assignments of original components.}
#' \item{n_iter}{Number of iterations until convergence.}
#' \item{converged}{Logical indicating whether convergence was achieved.}
#' \item{zeta_used}{The value of \code{zeta} used (either supplied or optimized).}
#'
#' @examples
#' orig_weights <- c(0.3, 0.3, 0.3, 0.1)
#' orig_alphas <- c(3, 4, 3.5, 7)
#' orig_betas  <- c(5, 6, 5.5, 2)
#' result <- BMR.CTD(orig_weights, orig_alphas, orig_betas, M = 2)
#' curve({dmix.beta(x, orig_weights, orig_alphas, orig_betas)}, from = 0, to = 1,
#'       col = "blue", lwd = 2, xlab = "x", ylab = "Density", main = "Original vs Reduced Mixture")
#' curve({dmix.beta(x, result$reduced_weights, result$reduced_alphas, result$reduced_betas)},
#'       add = TRUE, col = "red", lwd = 2, lty = 2)
#' @export
BMR.CTD <- function(orig_weights, orig_alphas, orig_betas, M, zeta = NULL, 
                    max_iter  = 100, tol = 1e-6, n_sample  = 10000) {
  
  compute_cost_matrix <- function(ra, rb, z) {
    hell <- beta.Hellinger(orig_alphas, orig_betas, 
                           matrix(rep(ra, each = length(orig_alphas)), ncol = M),
                           matrix(rep(rb, each = length(orig_alphas)), ncol = M))
    js   <- beta.JS(orig_alphas, orig_betas,
                    matrix(rep(ra, each = length(orig_alphas)), ncol = M),
                    matrix(rep(rb, each = length(orig_alphas)), ncol = M))
    z * hell + (1 - z) * js
  }
  
  reduction_phase <- function(z) {
    init   <- init_reduced_beta_mixture(orig_weights, orig_alphas, orig_betas, M, n_sample)
    p_old  <- c(init$alphas, init$betas, init$weights)
    
    for (i in seq_len(max_iter)) {
      ra <- p_old[1:M]
      rb <- p_old[(M + 1):(2 * M)]
      
      cost_matrix <- compute_cost_matrix(ra, rb, z)
      assignment  <- max.col(-cost_matrix)
      
      new_rw <- tabulate(assignment, nbins = M)
      new_rw <- as.numeric(sapply(seq_len(M), function(m) sum(orig_weights[assignment == m])))
      
      ab_mat <- beta.barycenter(assignment, M, orig_weights, orig_alphas, orig_betas, z)
      new_ra <- ab_mat[, 1]
      new_rb <- ab_mat[, 2]
      
      if (anyNA(new_ra)) new_ra[is.na(new_ra)] <- ra[is.na(new_ra)]
      if (anyNA(new_rb)) new_rb[is.na(new_rb)] <- rb[is.na(new_rb)]
      
      p_new <- c(new_ra, new_rb, new_rw)
      if (anyNA(p_new) || max(abs(p_new - p_old)) < tol) break
      
      p_old <- p_new
    }
    
    list(par = p_old, iter = i, converged = (i < max_iter))
  }
  
  if (is.null(zeta)) {
    zeta <- optimize(function(z) {
      res <- reduction_phase(z)
      ra  <- res$par[1:M]
      rb  <- res$par[(M + 1):(2 * M)]
      cost_matrix <- compute_cost_matrix(ra, rb, z)
      sum(orig_weights * apply(cost_matrix, 1, min))
    }, interval = c(0, 1))$minimum
  }
  
  res <- reduction_phase(zeta)
  pf  <- res$par
  
  ra <- pf[1:M]
  rb <- pf[(M + 1):(2 * M)]
  rw <- pf[(2 * M + 1):(3 * M)]
  
  cost_matrix <- compute_cost_matrix(ra, rb, zeta)
  assignment  <- max.col(-cost_matrix)
  total_cost  <- sum(orig_weights * apply(cost_matrix, 1, min))
  
  list(
    reduced_alphas  = ra,
    reduced_betas   = rb,
    reduced_weights = rw,
    total_cost      = total_cost,
    assignments     = assignment,
    n_iter          = res$iter,
    converged       = res$converged,
    zeta_used       = zeta
  )
}
