#' earthquake
#'
#' @description This data set contains the number of major earthquakes (magnitude 7 or greater) in the world from
#' 1900 to 2006. The data are available in Table 1.1 of Zucchini & MacDonald (2009). Zucchini &
#' MacDonald (2009) suggested using a Poisson mixture model to fit the data. The data frame has 107
#' rows and 1 column.
#'
#' @usage data(earthquake)
#' @format This data frame contains one column:
#' number: number of major earthquakes in a year.
#'
#' @references Zucchini W, MacDonald IL (2009). Hidden Markov models for time series: an introduction using
#' R. CRC Press.
#' 
#' @export
#' 
#' @examples data(earthquake)
#' out <- pmle.pois(unlist(earthquake),2,1)
#' plotmix.pois(unlist(earthquake), alpha= out[[1]], mu = out[[2]], extra.height=1.6)
"earthquake"