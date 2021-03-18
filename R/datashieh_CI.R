#' Function to compute CI around Shieh's effect size estimators
#'
#' @param Group.1 a (non-empty) numeric vector of data values.
#' @param Group.2 a (non-empty) numeric vector of data values.
#' @param conf.level confidence level of the interval
#' @param unbiased a logical variable indicating whether to compute the biased or unbiased estimator.
#' If TRUE, unbiased estimator is computed (Hedges' g or Hedges' g'). Otherwise, bias estimator is computed (Cohen's d or Cohen's d').
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param na.rm set whether Missing Values should be excluded (na.rm = TRUE) or not (na.rm = FALSE) - defaults to TRUE.
#'
#' @export datashieh_CI
#'
#' @exportS3Method datashieh_CI default
#' @exportS3Method print datashieh_CI
#'
#' @keywords Cohen's effect sizes, confidence interval
#' @return Returns Cohen's estimators of effect size and (1-alpha)% confidence interval around it, standard error
#' @importFrom stats na.omit sd pt uniroot

datashieh_CI <- function(Group.1,Group.2,conf.level,unbiased, alternative,na.rm) UseMethod("datashieh_CI")

datashieh_CIEst <- function(Group.1,
                            Group.2,
                            conf.level=.95,
                            unbiased=TRUE,
                            alternative="two.sided",
                            na.rm=TRUE){

  if (na.rm == TRUE ) {
    Group.1 <- na.omit(Group.1)
    Group.2 <- na.omit(Group.2)
  } else {
    Group.1 <- Group.1
    Group.2 <- Group.2
  }

  if(inherits(Group.1,c("numeric","integer")) == FALSE |inherits(Group.2,c("numeric","integer")) == FALSE)
    stop("Data are neither numeric nor integer")

    n1 <- length(Group.1)
    n2 <- length(Group.2)
    N <- n1+n2
    m1 <- mean(Group.1)
    m2 <- mean(Group.2)
    sd1 <- sd(Group.1)
    sd2 <- sd(Group.2)

    q1 <- n1/N
    q2 <- n2/N
    shieh.d <- (m1-m2)/sqrt(sd1^2/q1+sd2^2/q2)
    df <- ((sd1^2/n1+sd2^2/n2)^2)/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))
    w_obs <- (m1-m2)/sqrt(sd1^2/n1+sd2^2/n2)

    if(unbiased==TRUE){
      corr <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
    } else {corr <- 1}

    ES <- shieh.d*corr

    if(alternative=="two.sided"){

      # lower limit = limit of lambda such as 1-pt(q=t_obs, df=df, ncp = lambda) = (1-conf.level)/2 = alpha/2
      f=function(lambda,rep) 1-pt(q=w_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      lambda.1 <- out$root
      delta.1 <- lambda.1/sqrt(N) # lambda = delta * sqrt(N)
                                  # <--> delta = lambda/sqrt(N)

      # upper limit = limit of lambda such as pt(q=t_obs, df=df, ncp = lambda) = (1-conf.level)/2 = alpha/2
      f=function(lambda,rep) pt(q=w_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      lambda.2 <- out$root
      delta.2 <- lambda.2/sqrt(N) # lambda = delta * sqrt(N)
                                  # <--> delta = lambda/sqrt(N)

      result <- c(delta.1*corr, delta.2*corr)

    } else if (alternative == "greater"){

      # lower limit = limit of lambda such as 1-pt(q=t_obs, df=df, ncp = lambda) = (1-conf.level) = alpha
      f=function(lambda,rep) 1-pt(q=w_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      lambda.1 <- out$root
      delta.1 <- lambda.1/sqrt(N)

      # upper limit = + Inf
      delta.2 <- +Inf # if our expectation is mu1 > mu2, then we expect that (mu1-mu2)> 0 and therefore
                      # we want to check only the lower limit of the CI

      result <- c(delta.1*corr, delta.2)

    } else if (alternative == "less"){

      # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
      # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))

      delta.1 <- -Inf # if our expectation is mu1 < mu2, then we expect that (mu1-mu2)< 0 and therefore
                      # we want to check only the upper limit of the CI

      # upper limit = limit of lambda such as pt(q=t_obs, df=df, ncp = lambda) = (1-conf.level) = alpha
      f=function(lambda,rep) pt(q=w_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      lambda.2 <- out$root
      delta.2 <- lambda.2/sqrt(N) # See explanation in two.sided CI

      result <- c(delta.1, delta.2*corr)


    }

  # print results
  meth <- "Confidence interval around the raw mean difference"

  # Return results in list()
  invisible(
    list(ES = ES,
         conf.level = conf.level,
         CI = result)
  )

}

# Adding a default method in defining a function called datashieh_CI.default
datashieh_CI.default <- function(
  Group.1,
  Group.2,
  conf.level=.95,
  unbiased=TRUE,
  alternative="two.sided",
  na.rm=TRUE){

  out <- datashieh_CIEst(Group.1,Group.2,conf.level,unbiased,alternative,na.rm)
  out$ES <- out$ES
  out$call <- match.call()
  out$CI <- out$CI
  out$conf.level <- out$conf.level

  class(out) <- "datashieh_CI"
  out
}

print.datashieh_CI <- function(x,...){
  cat("Call:\n")
  print(x$call)

  cat("\nEffect size estimate :\n")
  print(round(x$ES,3))

  cat(paste0("\n",x$conf.level*100," % confidence interval around effect size estimate:\n"))
  print(round(x$CI,3))

}

