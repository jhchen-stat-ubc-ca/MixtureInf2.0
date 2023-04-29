#' SLC
#'
#' @description This data set contains 190 SLC measurements studied by Roeder, K. (1994). Roeder (1994) analyzed this data and concluded that a three component normal mixture with equal variance is most
#' suitable. Chen et al. (2012) also analyzed this data and gave a slightly better fit by a two component
#' normal mixture with unequal variance. The data frame has 190 rows and 1 column.
#'
#' @usage data(SLC)
#' @format This data frame contains one column:
#' SLC: 190 slc measurements.
#'
#' @references Chen, J., Li, P. and Fu, Y. (2012). Inference on the order of a normal mixture. JASA. 107, 1096-1105.
#' Roeder, K. (1994), A Graphical Technique for Determining the Number of Components in a Mixture of Normals, Journal of the American Statistical Association, 89, 487-500.

#' 
#' @examples data(SLC)
#' a <- rbind(c(0.6,0.4),c(0.2,0.3),c(0.01,0.01))
#' out1 <- pmle.norm(unlist(SLC),2,1,init.val = a)
#' out2 <- emtest.norm(unlist(SLC),m0 = 2,init.val = a)
#' plotmix.norm(unlist(SLC),alpha=out1[[1]][1,], mu=out1[[1]][2,], sigma=out1[[1]][3,],m0=2)
#' plotmix.norm(unlist(SLC),alpha=out2[[1]][1,], mu=out2[[1]][2,], sigma=out2[[1]][3,],m0=2)
"SLC"