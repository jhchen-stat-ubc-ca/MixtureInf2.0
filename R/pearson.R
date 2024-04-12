#' pearson
#'
#' @description This data set contains the ratio of "forehead" breadth to body length for 1000 crabs sampled at
#' Naples by Professor W.F.R. Weldon. Pearson (1894) used a two component normal mixture model
#' to fit this data set. The data frame has 29 rows and 2 columns.
#'
#' @usage data(pearson)
#' @format This data frame contains the following columns:
#' ratio: the boundaries of grouping intervals.
#' freq: the frequencies of observation falling into each interval.
#'
#' @references Pearson K (1894). Contributions to the mathematical theory of evolution. Philosophical Transactions of the Royal Society of London. A, 185, 71-110.
#' 
#' @examples data(pearson)
#' out <- pmle.norm(unlist(pearson),2,1)
#' plotmix.norm(unlist(pearson),alpha=out[[1]][1,], mu=out[[1]][2,], sigma=out[[1]][3,],m0=2)
#' @export
"pearson"