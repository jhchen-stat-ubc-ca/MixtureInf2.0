#' residual2
#'
#' @description This data set is from Lindsay and Roeder (1992). It includes the number of boys in families of size
#' 12. The number of families is 6115. The data frame has 13 rows and 2 columns.
#'
#' @usage data(residual2)
#' @format This data frame contains the 2 columns,
#' count: number of boys in family.
#' freq: number of families with corresponding.
#'
#' @references Lindsay, B. G. and Roeder, K. (1992). Residual diagnostics for mixture models. Journal of the
#' American Statistical Association,87(419), 785-794.
#' 
#' @example data(residual2)
#' pmle.binom(residual2, size = 12, m0 = 2,lambda = 1)
"residual2"