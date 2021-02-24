#' Function to compute CI around Cohen's effect size estimators
#'
#' @param m1 the average score of the first group
#' @param m2 the average score of the second group
#' @param sd1 the standard deviation the first group
#' @param sd2 the standard deviation the second group
#' @param n1 the first sample size
#' @param n2 the second sample size
#' @param conf.level confidence level of the interval
#' @param var.equal a logical variable indicating whether to assume equality of population variances.
#' If TRUE the pooled variance is used to estimate the standard error (= Cohen's d or Hedges' g). Otherwise, the square root of the non pooled
#' average of both variance estimates is used to estimate the standard error (Cohen's d' or Hedges' g').
#' @param unbiased a logical variable indicating whether to compute the biased or unbiased estimator.
#' If TRUE, unbiased estimator is computed (Hedges' g or Hedges' g'). Otherwise, bias estimator is computed (Cohen's d or Cohen's d').
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param na.rm set whether Missing Values should be excluded (na.rm = TRUE) or not (na.rm = FALSE) - defaults to TRUE.
#'
#' @export cohen.CI
#'
#' @keywords Cohen's effect sizes, confidence interval
#' @return Returns Cohen's estimators of effect size and (1-alpha)% confidence interval around it, standard error
#' @importFrom stats na.omit sd pt uniroot

cohen.CI <- function(m1,m2,sd1,sd2,n1,n2,conf.level,var.equal,unbiased, alternative,na.rm) UseMethod("cohen.CI")

cohen.CIEst <- function(m1,m2,sd1,sd2,n1,n2,
                        conf.level=.95,
                        var.equal=FALSE,
                        unbiased=TRUE,
                        alternative="two.sided",
                        na.rm=TRUE){

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

  if(var.equal==TRUE){

    pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
    t_obs <- (m1-m2)/sqrt(pooled_sd*(1/n1+1/n2))
    df <- n1+n2-2
    cohen.d <- (m1-m2)/pooled_sd

    if(unbiased==TRUE){
      corr <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
    } else {corr <- 1}

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
    df <- ((sd1^2/n1+sd2^2/n2)^2)/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))
    w_obs <- (m1-m2)/sqrt(sd1^2/n1+sd2^2/n2)

    # useful parameters (related with the link between Welch's t-statistic and Cohen's d's)
    s_bal <- sqrt((sd1^2+sd2^2)/2)
    N <- n1+n2
    q1 <- n1/N
    q2 <- n2/N
    s_unbal <- sqrt((1-q1)*sd1^2+(1-q2)*sd2^2)
    nratio <- n1/n2

    if(unbiased==TRUE){
      corr <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
    } else {corr <- 1}

    if(alternative=="two.sided"){

      # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level)/2 = alpha/2
      # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))

      f=function(lambda,rep) 1-pt(q=w_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      lambda.1 <- out$root
      delta.1 <- lambda.1/(sqrt(N)*(s_bal*sqrt(nratio))/((nratio+1)*s_unbal))
      # ncp_welch <- cohen.delta_prime *((sigma_bal*sqrt(nratio))/((nratio+1)*sigma_unbal))
      # with sigma_bal <- sqrt((sigma1^2+sigma2^2)/2)
      # and sigma_unbal <- sqrt((1-q1)*sigma1^2+(1-q2)*sigma2^2)
      # q1 <- n1/N and q2 <- n2/N
      # <--> cohen.delta_prime <- ncp_welch/(sqrt(N)*(sigma_bal*sqrt(nratio))/((nratio+1)*sigma_unbal))
      # sigma1 and sigma2 are unknown, so we estimate them with sd1 and sd2.

      # upper limit = limit of lambda such as pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level)/2 = alpha/2
      # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))
      f=function(lambda,rep) pt(q=w_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      lambda.2 <- out$root
      delta.2 <- lambda.2/(sqrt(N)*(s_bal*sqrt(nratio))/((nratio+1)*s_unbal))

      result <- c(delta.1*corr, delta.2*corr)

    } else if (alternative == "greater"){

      # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
      # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))

      f=function(lambda,rep) 1-pt(q=w_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      lambda.1 <- out$root
      delta.1 <- lambda.1/(sqrt(N)*(s_bal*sqrt(nratio))/((nratio+1)*s_unbal))

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

      f=function(lambda,rep) pt(q=w_obs, df=df, ncp = lambda)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      lambda.2 <- out$root
      delta.2 <- lambda.2/(sqrt(N)*(s_bal*sqrt(nratio))/((nratio+1)*s_unbal))

      result <- c(delta.1, delta.2*corr)


    }

  }

  # print results
  meth <- "Confidence interval around the raw mean difference"

  # Return results in list()
  invisible(
    list(cohen.d = cohen.d,
         conf.level = conf.level,
         CI = result)
  )

}

# Adding a default method in defining a function called cohen.CI.default
cohen.CI.default <- function(m1,m2,sd1,sd2,
                             n1,n2,conf.level=.95,
                             var.equal=FALSE,
                             unbiased=TRUE,
                             alternative="two.sided",
                             na.rm=TRUE){

  out <- cohen.CIEst(m1,m2,sd1,sd2,n1,n2,conf.level,var.equal,unbiased,alternative,na.rm)
  out$cohen.d <- out$cohen.d
  out$call <- match.call()
  out$CI <- out$CI
  out$conf.level <- out$conf.level

  class(out) <- "cohen.CI"
  out
}

print.cohen.CI <- function(x,...){
  cat("Call:\n")
  print(x$call)

  cat("\nRaw means difference :\n")
  print(round(x$cohen.d,3))

  cat(paste0("\n",x$conf.level*100," % confidence interval around the raw means difference:\n"))
  print(round(x$CI,3))

}


