#' pollen
#'
#' @description This data set contains the forest pollen counts from Mosimann 1962.The pollen data
#'  frame has four columns and 76 rows.The rows each sum to 100; the values are counts of four different types of pollen. Each row corresponds to a different level in the core; 
#'  the levels are in sequence with the first row being most recent and the last row being the oldest.
#'
#' @usage data(pollen)
#' @format This data frame contains four column:
#' Pinus
#' Abies
#' Quercus
#' Alnus
#'
#' @references J. E. Mosimann 1962. "On the compound multinomial distribution, the multivariate B-distribution,
#' and correlations among proportions". Biometrika, volume 49, numbers 1 and 2, pp65-82.
#' 
#' @examples data(pollen)
#' out1 <- pmle.norm(as.matrix(pollen),3,1)
#' out2 <- emtest.norm(as.matrix(pollen),m0 = 3)
#' plotmix.norm(as.matrix(pollen), alpha = out1[[1]][1,], mu = out1[[1]][2,], sigma = out1[[1]][3,], m0 = 3)
#' plotmix.norm(as.matrix(pollen), alpha = out2[[1]][1,], mu = out2[[1]][2,], sigma = out2[[1]][3,], m0 = 3)
"pollen"