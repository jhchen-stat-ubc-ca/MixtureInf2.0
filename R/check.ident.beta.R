#' check.comp
#'
#' @description The sub function that checks whether a group of components can be merged.
#' @param indices The index numbers of the input parameters.
#' @param mix_porp A vector of mixing proportions.
#' @param alpha A vector of subpopulation alphas.
#' @param beta A vector of subpopulation betas.
#' @param tol The tolerance value for checking the difference between actual and expected weights.
#' @export
check.comp <- function(indices, mix_porp, alpha, beta, tol = 1e-6) {
  a <- alpha[indices]
  b <- beta[indices]
  w <- mix_porp[indices]
  
  ord <- order(a)
  a <- a[ord]
  b <- b[ord]
  w <- w[ord]
  
  m <- length(indices)
  if (m > 1 && all(abs(diff(a)) < tol) && all(abs(diff(b)) < tol)) {
    return(list(a0 = a[1],
                b0 = b[1],
                n  = 0,
                indices = indices,
                new_weight = sum(w)))
  }
  
  Tval <- a[1] + b[1]
  if (any(abs((a + b) - Tval) > tol)) return(NULL)
  if (any(abs(diff(a) - 1) > tol)) return(NULL)
  
  a0 <- a[1]
  s  <- m - 1
  b0 <- Tval - a0 - s
  
  expected_b <- (Tval - a0) - (0:s)
  if (any(abs(b - expected_b) > tol)) return(NULL)
  
  exp_w <- choose(s, 0:s) * (beta(a0 + 0:s, b0 + s - 0:s) / beta(a0, b0))
  W <- sum(w)
  exp_w <- exp_w * (W / sum(exp_w))
  
  if (all(abs(w - exp_w) < tol))
    return(list(a0 = a0, b0 = b0, s = s, indices = indices, new_weight = W))
  else
    return(NULL)
}

#' iterative.reduce
#'
#' @description Repeatedly searches for any mergeable group among the current components.
#' When a mergeable group is found, they are merged (keeping track of original indices) and the search restarts.
#' @param mix_porp A vector of mixing proportions.
#' @param alpha A vector of subpopulation alphas.
#' @param beta A vector of subpopulation betas.
#' @param tol The tolerance value for checking the difference between actual and expected weights.
#' @export
iterative.reduce <- function(mix_porp, alpha, beta, tol = 1e-6) {
  comps <- lapply(seq_along(mix_porp), function(i) {
    list(weight = mix_porp[i], alpha = alpha[i], beta = beta[i], indices = i)
  })
  
  changed <- TRUE
  while (changed && length(comps) > 1) {
    changed <- FALSE
    n <- length(comps)
    weights <- sapply(comps, function(x) x$weight)
    alphas  <- sapply(comps, function(x) x$alpha)
    betas   <- sapply(comps, function(x) x$beta)
    
    merged <- NULL
    merged_comb <- NULL
    for (r in n:2) {
      for (cmb in combn(n, r, simplify = FALSE)) {
        candidate <- check.comp(seq_along(cmb),
                                weights[cmb],
                                alphas[cmb],
                                betas[cmb],
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
#' @description The main function that checks the identifiability of a beta mixture.
#' @param mix_porp A vector of mixing proportions.
#' @param alpha A vector of subpopulation alphas.
#' @param beta A vector of subpopulation betas.
#' @param tol The tolerance value for checking the difference between actual and expected weights.
#' @examples 
#' check.ident.beta(c(.03, .4, .04, .03, .25, .25), c(4, 5, 3, 2, 5.5, 6.5), c(2, 5, 3, 4, 6.5, 5.5))
#' @export
check.ident.beta <- function(mix_porp, alpha, beta, tol = 1e-6) {
  final_comps <- iterative.reduce(mix_porp, alpha, beta, tol)
  final_count <- length(final_comps)
  
  if (final_count == length(mix_porp))
    message("This beta mixture cannot be reduced.")
  else if (final_count == 1)
    message("This beta mixture distribution can be reduced to a beta distribution.")
  else
    message(sprintf("This beta mixture can be reduced to a %d component beta mixture.", final_count))
  
  comp_df <- do.call(rbind, lapply(seq_along(final_comps), function(i) {
    comp <- final_comps[[i]]
    data.frame(Component = i,
               Weight    = comp$weight,
               Alpha     = comp$alpha,
               Beta      = comp$beta,
               Combined_original_components = paste(comp$indices, collapse = ", "),
               stringsAsFactors = FALSE)
  }))
  
  return(comp_df)
}

