#' grains
#'
#' @description This data set contains the square root of the total number of grains for each planty from Loisel et al.,
#' (1994). Loisel et al., (1994) suggested using a finite normal mixture model. The grains data frame
#' has 150 rows and 1 column.
#'
#' @usage data(grains)
#' @format This data frame contains one column:
#' x: square root of the total number of grains for each planty.
#'
#' @references Loisel, P., Goffinet, B., Monod, H., and Montes De Oca, G. (1994). Detecting a major gene in an
#' F2 population. Biometrics, 50, 512–516.
#' 
#' @example data(grains)
#' pmle.norm(unlist(grains),2,1)
"grains"