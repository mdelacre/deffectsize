% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datashieh_CI.R
\name{datashieh_CI}
\alias{datashieh_CI}
\title{Function to compute CI around Shieh's effect size estimators}
\usage{
datashieh_CI(Group.1, Group.2, conf.level, unbiased, alternative, na.rm)
}
\arguments{
\item{Group.1}{a (non-empty) numeric vector of data values.}

\item{Group.2}{a (non-empty) numeric vector of data values.}

\item{conf.level}{confidence level of the interval}

\item{unbiased}{a logical variable indicating whether to compute the biased or unbiased estimator.
If TRUE, unbiased estimator is computed (Hedges' g or Hedges' g'). Otherwise, bias estimator is computed (Cohen's d or Cohen's d').}

\item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".}

\item{na.rm}{set whether Missing Values should be excluded (na.rm = TRUE) or not (na.rm = FALSE) - defaults to TRUE.}
}
\value{
Returns Cohen's estimators of effect size and (1-alpha)\% confidence interval around it, standard error
}
\description{
Function to compute CI around Shieh's effect size estimators
}
\keyword{Cohen's}
\keyword{confidence}
\keyword{effect}
\keyword{interval}
\keyword{sizes,}
