#' timesoffailure
#'
#' @description This data set contains the times of successive failures for the air conditioning system of each member
#' in a fleet of 13 Boeing 720 jet aircrafts. The data frame has 213 rows and 1 column.
#'
#' @usage data(timesoffailure)
#' @format This data frame contains one column,
#' x: pooled failure times about 213 observations.
#'
#' @references Proschan, F. (1963). Theoretical explanation of observed decreasing failure rate. Technometrics 5,
#' 375-83.
#' 
#' @examples data(timesoffailure)
#' out1 <- pmle.exp(unlist(timesoffailure),2,1)
#' out2 <- emtest.exp(unlist(timesoffailure),2)
#' plotmix.exp(alpha = out1[[1]], mu = out1[[2]], qq = 0.995, unlist(timesoffailure))
#' plotmix.exp(alpha = out2[[1]][1,], mu = out2[[1]][2,], qq = 0.995, unlist(timesoffailure))
#' @source <https://github.com/jhchen-stat-ubc-ca/Mixturelnf2.0>
"timesoffailure"