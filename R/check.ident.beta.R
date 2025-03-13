#' check_candidate
#'
#' @description # The sub function that actually checks whether a group of components can be merged or not. 
#' It is used in the check.ident.beta function.
#' @param indices The index number of the original input parameters.
#' @param mix_porp A vector of mixing proportions.
#' @param alpha A vector of subpopulation alphas
#' @param beta A vector of subpopulation betas
#' @param tol The tolerance value for checking the difference between actual and expected weights.
#' 
#' @export
check_candidate <- function(indices, mix_porp, alpha, beta, tol = 1e-6) {
  candidate_alphas  <- alpha[indices]
  candidate_betas   <- beta[indices]
  candidate_weights <- mix_porp[indices]
  
  o <- order(candidate_alphas)
  candidate_alphas  <- candidate_alphas[o]
  candidate_betas   <- candidate_betas[o]
  candidate_weights <- candidate_weights[o]
  
  m <- length(indices)
  
  if (m > 1 && all(abs(diff(candidate_alphas)) < tol) && all(abs(diff(candidate_betas)) < tol)) {
    return(list(a0 = candidate_alphas[1],
                b0 = candidate_betas[1],
                n = 0,
                indices = indices,
                new_weight = sum(candidate_weights)))
  }
  
  total_shapes <- candidate_alphas + candidate_betas
  Tval <- total_shapes[1]
  if (any(abs(total_shapes - Tval) > tol)) return(NULL)
  
  if (any(abs(diff(candidate_alphas) - 1) > tol)) return(NULL)
  
  a0 <- candidate_alphas[1]
  n  <- m - 1
  b0 <- Tval - a0 - n
  
  expected_betas <- (Tval - a0) - (0:n)
  if (any(abs(candidate_betas - expected_betas) > tol)) return(NULL)
  
  expected_weights <- sapply(0:n, function(k) {
    choose(n, k) * (beta(a0 + k, b0 + n - k) / beta(a0, b0))
  })
  
  W <- sum(candidate_weights)
  expected_weights <- expected_weights * (W / sum(expected_weights))
  
  if (all(abs(candidate_weights - expected_weights) < tol))
    return(list(a0 = a0, b0 = b0, n = n, indices = indices, new_weight = W))
  else
    return(NULL)
}
#' iterative_reduce
#'
#' @description # The sub function that repeatedly search for any mergeable group among the current components.
#' When a mergeable group is found, we merge them (keeping track of original indices) and restart.
#' It is used in the check.ident.beta function.
#' @param mix_porp A vector of mixing proportions.
#' @param alpha A vector of subpopulation alphas
#' @param beta A vector of subpopulation betas
#' @param tol The tolerance value for checking the difference between actual and expected weights.
#' 
#' @export
iterative_reduce <- function(mix_porp, alpha, beta, tol = 1e-6) {
  comps <- lapply(seq_along(mix_porp), function(i) {
    list(weight = mix_porp[i], alpha = alpha[i], beta = beta[i], indices = i)
  })
  
  changed <- TRUE
  while (changed && length(comps) > 1) {
    changed <- FALSE
    n <- length(comps)
    comp_weights <- sapply(comps, function(x) x$weight)
    comp_alpha   <- sapply(comps, function(x) x$alpha)
    comp_beta    <- sapply(comps, function(x) x$beta)
    
    merged <- NULL
    merged_comb <- NULL  
    
    for (r in n:2) {
      comb_list <- combn(n, r, simplify = FALSE)
      for (cmb in comb_list) {
        candidate <- check_candidate(seq_along(cmb),
                                     comp_weights[cmb],
                                     comp_alpha[cmb],
                                     comp_beta[cmb],
                                     tol)
        if (!is.null(candidate)) {
          merged <- candidate
          merged_comb <- cmb  
          break
        }
      }
      if (!is.null(merged)) break
    }
    
    if (!is.null(merged)) {
      new_indices <- sort(unlist(lapply(merged_comb, function(j) comps[[j]]$indices)))
      new_comp <- list(weight = merged$new_weight,
                       alpha  = merged$a0,
                       beta   = merged$b0,
                       indices = new_indices)
      comps <- comps[-merged_comb]
      comps <- c(comps, list(new_comp))
      changed <- TRUE
    }
  }
  return(comps)
}

#' check.ident.beta
#'
#' @description # The main function that checks the identifiability of a beta mixture.
#' @param mix_porp A vector of mixing proportions.
#' @param alpha A vector of subpopulation alphas
#' @param beta A vector of subpopulation betas
#' @param tol The tolerance value for checking the difference between actual and expected weights.
#' 
#' @examples 
#' check.ident.beta(c(.03,.4,.04,.03,.25,.25),c(4,5,3,2,5.5,6.5),c(2,5,3,4,6.5,5.5))
#' 
#' @export
check.ident.beta <- function(mix_porp, alpha, beta, tol = 1e-6) {
  final_comps <- iterative_reduce(mix_porp, alpha, beta, tol)
  final_count <- length(final_comps)
  
  if (final_count == length(mix_porp)) {
    message("This beta mixture cannot be reduced.")
  } else if (final_count == 1) {
    message("This beta mixture distribution can be reduced to a beta distribution.")
  } else {
    message(sprintf("This beta mixture can be reduced to a %d component beta mixture.", final_count))
  }
  
  comp_df <- do.call(rbind, lapply(seq_along(final_comps), function(i) {
    comp <- final_comps[[i]]
    data.frame(Component = i,
               Weight    = comp$weight,
               Alpha     = comp$alpha,
               Beta      = comp$beta,
               Combined_original_components   = paste(comp$indices, collapse = ", "),
               stringsAsFactors = FALSE)
  }))
  
  return(comp_df)
}