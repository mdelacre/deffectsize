% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datameandiff_CI.R
\name{datameandiff_CI}
\alias{datameandiff_CI}
\title{Function to compute CI around the raw mean difference}
\usage{
datameandiff_CI(Group.1, Group.2, conf.level, var.equal, alternative, na.rm)
}
\arguments{
\item{Group.1}{a (non-empty) numeric vector of data values.}

\item{Group.2}{a (non-empty) numeric vector of data values.}

\item{conf.level}{confidence level of the interval}

\item{var.equal}{a logical variable indicating whether to assume equality of population variances.
If TRUE the pooled variance is used to estimate the standard error. Otherwise, the standard error is estimated based on
unpooled variance.}

\item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".}

\item{na.rm}{set whether Missing Values should be excluded (na.rm = TRUE) or not (na.rm = FALSE) - defaults to TRUE.}
}
\value{
Returns raw mean difference, (1-alpha)\% confidence interval around mean difference, standard error
}
\description{
Function to compute CI around the raw mean difference
}
\keyword{confidence}
\keyword{difference,}
\keyword{interval}
\keyword{mean}
