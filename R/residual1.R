#' residual1
#'
#' @description This data set is from Lindsay and Roeder (1992). It includes the number of boys in families of size
#' 8. The number of families is 53680. The data frame has 9 rows and 2 columns.
#'
#' @usage data(residual1)
#' @format This data frame contains the 2 columns,
#' count: number of boys in family.
#' freq: number of families with corresponding.
#'
#' @references Lindsay, B. G. and Roeder, K. (1992). Residual diagnostics for mixture models. Journal of the
#' American Statistical Association,87(419), 785-794.
#' 
#' @examples data(residual1)
#' out <- pmle.binom(as.matrix(residual1), size = 8, m0 = 2,lambda = 1)
#' plotmix.binom(as.matrix(residual1), size = 8, alpha = out[[1]], theta = out[[2]])
#' @source <https://github.com/jhchen-stat-ubc-ca/Mixturelnf2.0>
"residual1"