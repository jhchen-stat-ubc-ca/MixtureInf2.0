#' lengthofstay
#'
#' @description This data set contains the length of stay (in days) for 469 geriatric patients in a psychiatric hospital in northeast London in 1991. Harrison and Millard (1991) used a mixture of two exponential
#' distributions in their analysis. The data frame has 469 rows and 1 column.
#'
#' @usage data(lengthofstay)
#' @format This data frame contains one column,
#' length: length of stay for 469 geriatric patients in a psychiatric hospital.
#'
#' @references Harrison G, Millard P (1991). Balancing acute and long-term care: the mathematics of throughput
#' in departments of geriatric medicine. Methods of information in medicine, 30(3), 221-228.
#' 
#' @example data(lengthofstay)
#' pmle.exp(unlist(lengthofstay),2,1)
"lengthofstay"