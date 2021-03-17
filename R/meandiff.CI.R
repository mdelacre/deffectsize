#' Function to compute CI around the raw mean difference
#'
#' @param m1 the average score of the first group
#' @param m2 the average score of the second group
#' @param sd1 the standard deviation the first group
#' @param sd2 the standard deviation the second group
#' @param n1 the first sample size
#' @param n2 the second sample size
#' @param conf.level confidence level of the interval
#' @param var.equal a logical variable indicating whether to assume equality of population variances.
#' If TRUE the pooled variance is used to estimate the standard error. Otherwise, the standard error is estimated based on
#' unpooled variance.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#'
#' @export meandiff.CI
#'
#' @keywords mean difference, confidence interval
#' @return Returns raw mean difference, (1-alpha)% confidence interval around mean difference, standard error
#' @importFrom stats na.omit sd pt uniroot

meandiff.CI <- function(m1,m2,sd1,sd2,n1,n2,conf.level,var.equal,alternative) UseMethod("meandiff.CI")

meandiff.CIEst <- function(m1,m2,sd1,sd2,
                            n1,n2,conf.level=.95,
                            var.equal=FALSE,
                           alternative="two.sided"){

#Alert message if at least one of the parameter is neither numeric nor integer

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

      theo_mudiff.2 <- Inf  # if our expectation is mu1 > mu2, then we expect that (mu1-mu2)> 0 and therefore
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

      theo_mudiff.2 <-  Inf   # if our expectation is mu1 > mu2, then we expect that (mu1-mu2)> 0 and therefore
                              # we want to check only the lower limit of the CI

      result <- c(theo_mudiff.1, theo_mudiff.2)
    }


  } else if (alternative == "less"){

    if (var.equal==TRUE){
      pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
      SE <- pooled_sd*sqrt(1/n1+1/n2)


      # Lower limit = - inf

      theo_mudiff.1 <- -Inf  # if our expectation is mu1 < mu2, then we expect that (mu1-mu2)< 0 and therefore
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

# Adding a default method in defining a function called meandiff.CI.default

meandiff.CI.default <- function(m1,m2,sd1,sd2,
                           n1,n2,conf.level=.95,
                           var.equal=FALSE,
                           alternative="two.sided"){

  out <- meandiff.CIEst(m1,m2,sd1,sd2,n1,n2,conf.level,var.equal,alternative)
  out$Meandiff <- out$Meandiff
  out$std.error <- out$std.error
  out$call <- match.call()
  out$CI <- out$CI
  out$conf.level <- out$conf.level

  class(out) <- "meandiff.CI"
  out
}

print.meandiff.CI <- function(x,...){
  cat("Call:\n")
  print(x$call)

  cat("\nRaw means difference :\n")
  print(round(x$Meandiff,3))

  cat(paste0("\n",x$conf.level*100," % confidence interval around the raw means difference:\n"))
  print(round(x$CI,3))

}




