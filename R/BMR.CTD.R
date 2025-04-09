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
  H2b <- 1 - beta((a1 + a2) / 2, (b1 + b2) / 2) / sqrt(beta(a1, b1) * beta(a2, b2))
  return(H2b)  
}

#' beta.barycenter
#'
#' @description Computes the barycenter of a collection of beta distributions 
#' with respect to the Hellinger distance, using a gradient-based optimization method.
#' It is a sub function used in the BMR.CTD function.
#'
#' @param cluster_weights A numeric vector of weights associated with each beta distribution.
#' @param cluster_alphas A numeric vector of shape parameters \eqn{\alpha} for each beta distribution.
#' @param cluster_betas A numeric vector of shape parameters \eqn{\beta} for each beta distribution.
#'
#' @return A numeric vector of length two containing the estimated \code{(alpha, beta)} parameters 
#'         of the barycenter.
#' @export


beta.barycenter <- function(cluster_weights, cluster_alphas, cluster_betas) {
  eval_f_grad <- function(par) {
    a <- par[1]
    b <- par[2]
    
    obj <- sum(cluster_weights * beta.Hellinger(cluster_alphas, cluster_betas, a, b))
    
    u <- (cluster_alphas + a) / 2
    v <- (cluster_betas  + b) / 2
    
    dB_da <- 0.5 * beta(u, v) * (digamma(u) - digamma(u + v))
    dB_db <- 0.5 * beta(u, v) * (digamma(v) - digamma(u + v))
    
    dG_da <- 0.5 * sqrt(beta(a, b)) * (digamma(a) - digamma(a + b))
    dG_db <- 0.5 * sqrt(beta(a, b)) * (digamma(b) - digamma(a + b))
    
    
    dT_da <- (dB_da * sqrt(beta(a, b)) - beta(u, v) * dG_da) / (beta(a, b))
    dT_db <- (dB_db * sqrt(beta(a, b)) - beta(u, v) * dG_db) / (beta(a, b))
    
    grad_a <- -sum(cluster_weights * dT_da / sqrt(beta(cluster_alphas, cluster_betas)))
    grad_b <- -sum(cluster_weights * dT_db / sqrt(beta(cluster_alphas, cluster_betas)))
    
    return(list(objective = obj, gradient = c(grad_a, grad_b)))
  }
  
  w_mean <- sum(cluster_weights * (cluster_alphas / (cluster_alphas + cluster_betas)))
  a0 <- w_mean 
  b0 <- 1 - w_mean
  init_par <- c(a0, b0)
  
  res <- optim(
    par = init_par,
    fn = function(x) eval_f_grad(x)$objective,
    gr = function(x) eval_f_grad(x)$gradient,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    control = list(maxit = 1000, factr = 1e6)
  )
  
  return(c(res$par[1], res$par[2]))
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
  
  init_alphas <- numeric(M)
  init_betas  <- numeric(M)
  init_weights <- numeric(M)
  
  for(m in 1:M) {
    cluster_samples <- samples[km$cluster == m]
    cluster_mean <- mean(cluster_samples)
    cluster_var  <- var(cluster_samples)
    s <- cluster_mean * (1 - cluster_mean) / cluster_var - 1
    if(s <= 0) s <- 5
    init_alphas[m] <- cluster_mean * s
    init_betas[m] <- (1 - cluster_mean) * s
    init_weights[m] <- length(cluster_samples) / n_sample
  }
  return(list(alphas = init_alphas, betas = init_betas, weights = init_weights))
}


#' BMR.CTD
#'
#' @description performs beta mixture reduction by iteratively assigning components based on Hellinger distance 
#' and recomputing cluster barycenters.
#'
#' @param orig_weights A numeric vector of original mixture weights.
#' @param orig_alphas A numeric vector of \eqn{\alpha} parameters for the original beta components.
#' @param orig_betas A numeric vector of \eqn{\beta} parameters for the original beta components.
#' @param M The number of components in the reduced mixture.
#' @param max_iter Maximum number of iterations for the algorithm. Default is 100.
#' @param tol Tolerance for convergence based on change in total cost. Default is 1e-6.
#' @param n_sample Number of samples to use for initialization. Default is 10000.
#'
#' @return A list with elements:
#' \item{reduced_alphas}{Estimated \eqn{\alpha} parameters of the reduced mixture.}
#' \item{reduced_betas}{Estimated \eqn{\beta} parameters of the reduced mixture.}
#' \item{reduced_weights}{Weights of the reduced mixture components.}
#' \item{cost_history}{Vector of objective costs at each iteration.}
#' \item{assignments}{Cluster assignment for each original component.}
#' \item{n_iter}{Number of iterations run.}
#' \item{converged}{Logical indicating whether convergence was achieved.}
#' 
#' @examples
#' orig_weights <- c(0.3, 0.3, 0.3, 0.1)
#' orig_alphas <- c(3, 4, 3.5, 7)
#' orig_betas  <- c(5, 6, 5.5, 2)
#' result <- BMR.CTD(orig_weights, orig_alphas, orig_betas, M = 2,
#' max_iter = 100, tol = 1e-6, n_sample = 10000)
#' curve({dmix.beta(x,orig_weights,orig_alphas,orig_betas)}, from = 0, to = 1, col = "blue", lwd = 2,
#' xlab = "x", ylab = "Density", main = "Original vs Reduced Beta Mixture")
#' 
#' curve({dmix.beta(x,result$reduced_weights,result$reduced_alphas,result$reduced_betas)}, 
#' add = TRUE, col = "red", lwd = 2, lty = 2)
#' @export

BMR.CTD <- function(orig_weights, orig_alphas, orig_betas, M, 
                    max_iter = 100, tol = 1e-6, n_sample = 10000) {
  N <- length(orig_weights)
  
  init_res <- init_reduced_beta_mixture(orig_weights, orig_alphas, orig_betas, M, n_sample)
  red_alphas <- init_res$alphas
  red_betas  <- init_res$betas
  red_weights <- init_res$weights
  
  cost_history <- numeric(max_iter)  
  converged <- FALSE
  assignment <- rep(NA_integer_, N)
  
  for(iter in 1:max_iter) {
    # Majorization step
    cost_matrix <- sapply(1:M, function(m) {
      beta.Hellinger(orig_alphas, orig_betas, red_alphas[m], red_betas[m])
    })
    assignment <- max.col(-cost_matrix)  
    cost_vec <- apply(cost_matrix, 1, min)
    total_cost <- sum(orig_weights * cost_vec)
    cost_history[iter] <- total_cost
    
    new_red_weights <- sapply(1:M, function(m) sum(orig_weights[assignment == m]))
    
    new_red_alphas <- red_alphas
    new_red_betas  <- red_betas
    for(m in 1:M) {
      idx <- which(assignment == m)
      if(length(idx) > 0) {
        # Minimization step
        new_params <- beta.barycenter(orig_weights[idx], orig_alphas[idx], orig_betas[idx])
        new_red_alphas[m] <- new_params[1]
        new_red_betas[m]  <- new_params[2]
      }
    }
    
    if(iter > 1 && abs(cost_history[iter-1] - total_cost) < tol) {
      converged <- TRUE
      cost_history <- cost_history[1:iter]
      break
    }
    
    red_alphas <- new_red_alphas
    red_betas  <- new_red_betas
    red_weights <- new_red_weights
  }
  
  return(list(reduced_alphas = red_alphas,
              reduced_betas  = red_betas,
              reduced_weights = red_weights,
              cost_history   = cost_history,
              assignments    = assignment,
              n_iter         = iter,
              converged      = converged))
}
