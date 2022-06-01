#' Function to compute CI around shieh's effect size estimators
#'
#' @param m1 the average score of the first group
#' @param m2 the average score of the second group
#' @param sd1 the standard deviation the first group
#' @param sd2 the standard deviation the second group
#' @param n1 the first sample size
#' @param n2 the second sample size
#' @param conf.level confidence level of the interval
#' @param unbiased a logical variable indicating whether to compute the biased or unbiased estimator.
#' If TRUE, unbiased estimator is computed (Hedges' g or Hedges' g'). Otherwise, bias estimator is computed (Cohen's d or Cohen's d').
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#'
#' @export shieh_CI
#'
#' @exportS3Method shieh_CI default
#' @exportS3Method print shieh_CI
#'
#' @keywords Cohen's effect sizes, confidence interval
#' @return Returns Cohen's estimators of effect size and (1-alpha)% confidence interval around it, standard error
#' @importFrom stats na.omit sd pt uniroot

shieh_CI <- function(m1,m2,sd1,sd2,n1,n2,conf.level,unbiased, alternative) UseMethod("shieh_CI")

shieh_CIEst <- function(m1,m2,sd1,sd2,n1,n2,
                            conf.level=.95,
                            unbiased=TRUE,
                            alternative="two.sided"){

  param <- data.frame(m1,m2,sd1,sd2,n1,n2)
  vect <- NULL
  for (i in seq_len(length(param))){
    if(inherits(param[,i],c("numeric","integer"))==FALSE){
      vect <-   c(vect,names(param[i]))
    } else {vect=vect}
  }

  if(inherits(c(m1,m2,sd1,sd2,n1,n2),c("numeric","integer"))==FALSE){
    if(length(vect)==1){
      obj <- vect
      alert="is neither numeric nor integer"
    } else if (length(vect)>1){
      obj <-paste(paste(vect[-length(vect)],collapse=", "),"and",vect[length(vect)])
      alert="are neither numeric nor integer"
    }

    stop(paste(obj,alert))

  }

  N <- n1+n2
  q1 <- n1/N
  q2 <- n2/N
  shieh.d <- (m1-m2)/sqrt(sd1^2/q1+sd2^2/q2)
  df <- ((sd1^2/n1+sd2^2/n2)^2)/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))
  w_obs <- (m1-m2)/sqrt(sd1^2/n1+sd2^2/n2)

  if(unbiased==TRUE){
    corr <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
  } else {corr <- 1}

  if(corr=="NaN"){
    alert2="Correction for bias is only for small sample sizes. Use 'unbiased=FALSE'"
    stop(alert2)
  } else {ES <- shieh.d*corr}

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
    delta.1 <- lambda.1/sqrt(N) # See explanation in two.sided CI

    # upper limit = + Inf
    delta.2 <- +Inf # if our expectation is mu1 > mu2, then we expect that (mu1-mu2)> 0 and therefore
                    # we want to check only the lower limit of the CI

    result <- c(delta.1*corr, delta.2)

  } else if (alternative == "less"){

    # lower limit = -Inf
    delta.1 <- -Inf  # if our expectation is mu1 < mu2, then we expect that (mu1-mu2)< 0 and therefore
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

# Adding a default method in defining a function called shieh_CI.default
shieh_CI.default <- function(
  m1,m2,sd1,sd2,n1,n2,
  conf.level=.95,
  unbiased=TRUE,
  alternative="two.sided"){

  out <- shieh_CIEst(m1,m2,sd1,sd2,n1,n2,
                     conf.level,unbiased,alternative)
  out$ES <- out$ES
  out$call <- match.call()
  out$CI <- out$CI
  out$conf.level <- out$conf.level

  class(out) <- "shieh_CI"
  out
}

print.shieh_CI <- function(x,...){
  cat("Call:\n")
  print(x$call)

  cat("\nEffect size estimate :\n")
  print(round(x$ES,3))

  cat(paste0("\n",x$conf.level*100," % confidence interval around effect size estimate:\n"))
  print(round(x$CI,3))

}
