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
#' @example data(timesoffailure)
#' pmle.exp(unlist(timesoffailure),2,1)
"timesoffailure"