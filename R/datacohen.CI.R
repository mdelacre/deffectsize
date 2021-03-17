#' Function to compute CI around Cohen's effect size estimators
#'
#' @param Group.1 a (non-empty) numeric vector of data values.
#' @param Group.2 a (non-empty) numeric vector of data values.
#' @param conf.level confidence level of the interval
#' @param var.equal a logical variable indicating whether to assume equality of population variances.
#' If TRUE the pooled variance is used to estimate the standard error (= Cohen's d or Hedges' g). Otherwise, the square root of the non pooled
#' average of both variance estimates is used to estimate the standard error (Cohen's d' or Hedges' g').
#' @param unbiased a logical variable indicating whether to compute the biased or unbiased estimator.
#' If TRUE, unbiased estimator is computed (Hedges' g or Hedges' g'). Otherwise, bias estimator is computed (Cohen's d or Cohen's d').
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param na.rm set whether Missing Values should be excluded (na.rm = TRUE) or not (na.rm = FALSE) - defaults to TRUE.
#'
#' @export datacohen.CI
#'
#' @keywords Cohen's effect sizes, confidence interval
#' @return Returns Cohen's estimators of effect size and (1-alpha)% confidence interval around it, standard error
#' @importFrom stats na.omit sd pt uniroot

datacohen.CI <- function(Group.1,Group.2,conf.level,var.equal,unbiased, alternative,na.rm) UseMethod("datacohen.CI")

datacohen.CIEst <- function(Group.1,
                               Group.2,
                               conf.level=.95,
                               var.equal=FALSE,
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

  if(var.equal==TRUE){

    pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
    t_obs <- (m1-m2)/sqrt(pooled_sd^2*(1/n1+1/n2))
    df <- n1+n2-2
    cohen.d <- (m1-m2)/pooled_sd

    if(unbiased==TRUE){
      corr <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
    } else {corr <- 1}

    ES <- cohen.d*corr

    if(alternative=="two.sided"){

      # lower limit = limit of lambda such as 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level)/2 = alpha/2
      f=function(lambda,rep) 1-pt(q=t_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      lambda.1 <- out$root
      delta.1 <- lambda.1*sqrt(1/n1+1/n2)   # See explanation for delta.2

      # upper limit = limit of lambda such as pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level)/2 = alpha/2
      f=function(lambda,rep) pt(q=t_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      lambda.2 <- out$root
      delta.2 <- lambda.2*sqrt(1/n1+1/n2)   # Pr[t_obs(alpha/2,df=df,ncp=lambda) <= lambda <= t_obs(1-alpha/2,df=df,ncp=lambda)=.025
      # because lambda = delta * sqrt[n1n2/(n1+n2)] :
      # Pr[t_obs(alpha/2,df=df,ncp=lambda) <= delta * sqrt[n1n2/(n1+n2)] <= t_obs(1-alpha/2,df=df,ncp=lambda)=.025
      # Pr[t_obs(alpha/2,df=df,ncp=lambda)*sqrt[(n1+n2)/n1n2] <= delta <= t_obs(1-alpha/2,df=df,ncp=lambda)*sqrt[(n1+n2)/n1n2]=.025
      # and sqrt[(n1+n2)/n1n2]=sqrt(1/n1+1/n2)
      result <- c(delta.1*corr, delta.2*corr)

    } else if (alternative == "greater"){

      # lower limit = limit of lambda such as 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level) = alpha
      f=function(lambda,rep) 1-pt(q=t_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      lambda.1 <- out$root
      delta.1 <- lambda.1*sqrt(1/n1+1/n2)   # See explanation for delta.2

      # upper limit = +Inf
      delta.2 <- +Inf
      result <- c(delta.1*corr, delta.2)

    } else if (alternative == "less"){

      # lower limit = -Inf
      delta.1 <- -Inf

      # upper limit = limit of lambda such as pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level) = alpha
      f=function(lambda,rep) pt(q=t_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      lambda.2 <- out$root
      delta.2 <- lambda.2*sqrt(1/n1+1/n2)

      result <- c(delta.1, delta.2*corr)

    }

  } else if (var.equal==FALSE){

    cohen.d <- (m1-m2)/sqrt((sd1^2+sd2^2)/2)
    df <- ((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
    t_obs <- (sqrt(n1*n2)*(m1-m2))/sqrt(n2*sd1^2+n1*sd2^2)

    if(unbiased==TRUE){
      corr <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
    } else {corr <- 1}

    ES <- cohen.d*corr
    if(alternative=="two.sided"){

      # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level)/2 = alpha/2
      f=function(lambda,rep) 1-pt(q=t_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      lambda.1 <- out$root
      delta.1 <- lambda.1*sqrt((2*(n2*sd1^2+n1*sd2^2))/(n1*n2*(sd1^2+sd2^2)))

      # upper limit = limit of lambda such as pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level)/2 = alpha/2
      f=function(lambda,rep) pt(q=t_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      lambda.2 <- out$root
      delta.2 <- lambda.2*sqrt((2*(n2*sd1^2+n1*sd2^2))/(n1*n2*(sd1^2+sd2^2)))

      result <- c(delta.1*corr, delta.2*corr)

    } else if (alternative == "greater"){

      # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
      # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))

      f=function(lambda,rep) 1-pt(q=t_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      lambda.1 <- out$root
      delta.1 <- lambda.1*sqrt((2*(n2*sd1^2+n1*sd2^2))/(n1*n2*(sd1^2+sd2^2)))

      # upper limit = limit of lambda such as pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
      # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))

      delta.2 <- +Inf

      result <- c(delta.1*corr, delta.2)

    } else if (alternative == "less"){

      # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
      # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))

      delta.1 <- -Inf

      # upper limit = limit of lambda such as pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
      # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))

      f=function(lambda,rep) pt(q=t_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      lambda.2 <- out$root
      delta.2 <- lambda.2*sqrt((2*(n2*sd1^2+n1*sd2^2))/(n1*n2*(sd1^2+sd2^2)))

      result <- c(delta.1, delta.2*corr)


    }

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

# Adding a default method in defining a function called datacohen.CI.default
datacohen.CI.default <- function(
  Group.1,
  Group.2,
  conf.level=.95,
  var.equal=FALSE,
  unbiased=TRUE,
  alternative="two.sided",
  na.rm=TRUE){

  out <- datacohen.CIEst(Group.1,Group.2,conf.level,var.equal,unbiased,alternative,na.rm)
  out$ES <- out$ES
  out$call <- match.call()
  out$CI <- out$CI
  out$conf.level <- out$conf.level

  class(out) <- "datacohen.CI"
  out
}

print.datacohen.CI <- function(x,...){
  cat("Call:\n")
  print(x$call)

  cat("\nEffect size estimate :\n")
  print(round(x$ES,3))

  cat(paste0("\n",x$conf.level*100," % confidence interval around effect size estimate:\n"))
  print(round(x$CI,3))

}

