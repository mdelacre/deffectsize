#' Function to compute CI around the raw mean difference
#'
#' @param Group.1 a (non-empty) numeric vector of data values.
#' @param Group.2 a (non-empty) numeric vector of data values.
#' @param conf.level confidence level of the interval
#' @param var.equal a logical variable indicating whether to assume equality of population variances.
#' If TRUE the pooled variance is used to estimate the standard error. Otherwise, the standard error is estimated based on
#' unpooled variance.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param na.rm set whether Missing Values should be excluded (na.rm = TRUE) or not (na.rm = FALSE) - defaults to TRUE.
#'
#' @export datameandiff.CI
#'
#' @exportS3Method datameandiff.CI default
#' @exportS3Method print datameandiff.CI
#'
#' @keywords mean difference, confidence interval
#' @return Returns raw mean difference, (1-alpha)% confidence interval around mean difference, standard error
#' @importFrom stats na.omit sd pt uniroot

datameandiff.CI <- function(Group.1,Group.2,conf.level,var.equal,alternative,na.rm) UseMethod("datameandiff.CI")

datameandiff.CIEst <- function(Group.1,
                            Group.2,
                            conf.level=.95,
                            var.equal=FALSE,
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
  m1 <- mean(Group.1)
  m2 <- mean(Group.2)
  sd1 <- sd(Group.1)
  sd2 <- sd(Group.2)

  if(alternative=="two.sided"){

    if (var.equal==TRUE){
      pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
      SE <- pooled_sd*sqrt(1/n1+1/n2) # standard error

      # lower limit = limit of mu1-mu2 such as 1-pt(q=t_obs, df=df) = (1-conf.level)/2 = alpha/2
      # with t_obs = ((m1-m2)-(theo_mudiff))/SE
      f=function(theo_mudiff,rep) 1-pt(q=((m1-m2)-(theo_mudiff))/SE, df=n1+n2-2)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      theo_mudiff.1 <- out$root

      # upper limit = limit of mu1-mu2 such aspt(q=t_obs, df=df) = (1-conf.level)/2 = alpha/2
      # with t_obs = ((m1-m2)-(theo_mudiff))/SE
      f=function(theo_mudiff,rep) pt(q=((m1-m2)-(theo_mudiff))/SE, df=n1+n2-2)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      theo_mudiff.2 <- out$root

      result <- c(theo_mudiff.1, theo_mudiff.2)

    } else if (var.equal==FALSE){

      SE <- sqrt(sd1^2/n1+sd2^2/n2) # standard error
      DF <- ((sd1^2/n1+sd2^2/n2)^2)/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))

      # lower limit = limit of mu1-mu2 such as 1-pt(q=t_obs, df=df) = (1-conf.level)/2 = alpha/2
      # with t_obs = ((m1-m2)-(theo_mudiff))/SE
      f=function(theo_mudiff,rep) 1-pt(q=((m1-m2)-(theo_mudiff))/SE, df=DF)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")

      theo_mudiff.1 <- out$root

      # upper limit = limit of mu1-mu2 such aspt(q=t_obs, df=df) = (1-conf.level)/2 = alpha/2
      # with t_obs = ((m1-m2)-(theo_mudiff))/SE
      f=function(theo_mudiff,rep) pt(q=((m1-m2)-(theo_mudiff))/SE, df=DF)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      theo_mudiff.2 <- out$root

      result <- c(theo_mudiff.1, theo_mudiff.2)
    }

  } else if (alternative == "greater"){

    if (var.equal==TRUE){
      pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
      SE <- pooled_sd*sqrt(1/n1+1/n2)


      # Lower limit = limit of mu1-mu2 such as 1-pt(q=t_obs, df=df) = (1-conf.level) = alpha
      # with t_obs = ((m1-m2)-(theo_mudiff))/SE
      f=function(theo_mudiff,rep) 1-pt(q=((m1-m2)-(theo_mudiff))/SE, df=n1+n2-2)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")

      theo_mudiff.1 <- out$root

      # upper limit = +Inf

      theo_mudiff.2 <- Inf # if our expectation is mu1 > mu2, then we expect that (mu1-mu2)> 0 and therefore
                           # we want to check only the lower limit of the CI

      result <- c(theo_mudiff.1, theo_mudiff.2)

    } else if (var.equal==FALSE) {
      SE <- sqrt(sd1^2/n1+sd2^2/n2)
      DF <- ((sd1^2/n1+sd2^2/n2)^2)/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))

      # Lower limit = limit of mu1-mu2 such as 1-pt(q=t_obs, df=DF) = (1-conf.level) = alpha
      # with t_obs = ((m1-m2)-(theo_mudiff))/SE
      f=function(theo_mudiff,rep) 1-pt(q=((m1-m2)-(theo_mudiff))/SE, df=DF)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      theo_mudiff.1 <- out$root

      # Upper limit = +Inf
      theo_mudiff.2 <-  Inf # if our expectation is mu1 > mu2, then we expect that (mu1-mu2)> 0 and therefore
                            # we want to check only the lower limit of the CI

      result <- c(theo_mudiff.1, theo_mudiff.2)
    }


  } else if (alternative == "less"){

    if (var.equal==TRUE){
      pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
      SE <- pooled_sd*sqrt(1/n1+1/n2)


      # Lower limit = - inf
      theo_mudiff.1 <- -Inf # if our expectation is mu1 < mu2, then we expect that (mu1-mu2)< 0 and therefore
                            # we want to check only  the upper limit of the CI

      # Upper limit = limit of mu1-mu2 such as pt(q=t_obs, df=df) = (1-conf.level) = alpha
      # with t_obs = ((m1-m2)-(theo_mudiff))/SE
      f=function(theo_mudiff,rep) pt(q=((m1-m2)-(theo_mudiff))/SE, df=n1+n2-2)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      theo_mudiff.2 <- out$root

      result <- c(theo_mudiff.1, theo_mudiff.2)

    } else if (var.equal==FALSE) {
      SE <- sqrt(sd1^2/n1+sd2^2/n2)
      DF <- ((sd1^2/n1+sd2^2/n2)^2)/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))

      # Lower limit = - inf
      theo_mudiff.1 <- -Inf # if our expectation is mu1 < mu2, then we expect that (mu1-mu2)< 0 and therefore
                            # we want to check only  the upper limit of the CI

      # Upper limit = limit of mu1-mu2 such as pt(q=t_obs, df=DF) = (1-conf.level) = alpha
      # with t_obs = ((m1-m2)-(theo_mudiff))/SE
      f=function(theo_mudiff,rep) pt(q=((m1-m2)-(theo_mudiff))/SE, df=DF)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      theo_mudiff.2 <- out$root

      result <- c(theo_mudiff.1, theo_mudiff.2)
    }

  }

  # print results
  meth <- "Confidence interval around the raw mean difference"

  # Return results in list()
  invisible(
    list(Meandiff = m1-m2,
         conf.level = conf.level,
         Std.error = SE,
         CI = result)
  )

}

# Adding a default method in defining a function called datameandiff.CI.default
datameandiff.CI.default <- function(
           Group.1,
           Group.2,
           conf.level=.95,
           var.equal=FALSE,
           alternative="two.sided",
           na.rm=TRUE){

  out <- datameandiff.CIEst(Group.1,Group.2,conf.level,var.equal,alternative,na.rm)
  out$Meandiff <- out$Meandiff
  out$std.error <- out$std.error
  out$call <- match.call()
  out$CI <- out$CI
  out$conf.level <- out$conf.level

  class(out) <- "datameandiff.CI"
  out
}

print.datameandiff.CI <- function(x,...){
  cat("Call:\n")
  print(x$call)

  cat("\nRaw means difference :\n")
  print(round(x$Meandiff,3))

  cat(paste0("\n",x$conf.level*100," % confidence interval around the raw means difference:\n"))
  print(round(x$CI,3))

}




