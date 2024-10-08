% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agesch.R
\docType{data}
\name{agesch}
\alias{agesch}
\title{agesch}
\format{
This data frame contains one column:
age: age of onset of schizophrenia of 152 males.
}
\source{
<https://github.com/jhchen-stat-ubc-ca/Mixturelnf2.0>
}
\usage{
data(agesch)
}
\description{
This data set is from a schizophrenia study by Levine (1981), who collected the results of seven 
studies on the age of onset of schizophrenia of 99 females and 152 males. This data set contains
observed values of only male data. Everitt et al., (2001) fitted the male data using a two-component
normal mixture model. The data frame has 152 rows and 1 column.
}
\examples{
data(agesch)
out <- pmle.norm(unlist(agesch),2,1)
plotmix.norm(unlist(agesch),alpha=out[[1]][1,], mu=out[[1]][2,], sigma=out[[1]][3,],m0=2)
}
\references{
Levine, R. (1981). Sex differences in schizophrenia: Timing or subtypes? Psychological Bulletin,
90, 432-444.
Everitt, B. S., Landau, S., and Leese, M. (2001). Cluster Analysis (4th ed.). New York, NY: Oxford
University Press, Inc.
}
\keyword{datasets}
