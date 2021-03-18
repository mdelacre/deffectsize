n1 <- 10
n2 <- 12
vd <- c(rnorm(n1,4,3),rnorm(n2,12,3))
vi <- c(rep(1,n1),rep(2,n2))
bdd <- data.frame(vi,vd)
Group.1 <- vd[vi==1]
Group.2 <- vd[vi==2]
m1 <- mean(Group.1)
m2 <- mean(Group.2)
sd1 <- sd(Group.1)
sd2 <- sd(Group.2)

res <- shieh(m1,m2,sd1,sd2,n1,n2,conf.level=.95,unbiased=F, alternative="two.sided")

test_that("data types correct", {
  expect_is(res,"shieh")
})

# Checking if both cohen and datacohen are consistent

Group.1 <- rnorm(10)
Group.2 <- rnorm(12)
conf.level=.90
alternative= "less"
na.rm=T
unbiased=F

m1 <-  mean(Group.1)
m2 <-  mean(Group.2)
n1 <- length(Group.1)
n2 <- length(Group.2)
sd1 <- sd(Group.1)
sd2 <- sd(Group.2)

res2 <- shieh(m1,m2,sd1,sd2,n1,n2,conf.level,unbiased, alternative)
res1 <- datashieh(Group.1,Group.2,conf.level,unbiased, alternative,na.rm)

testthat::test_that("functions are consistent",{
  expect_equal(res1$Meandiff,res2$Meandiff)
  expect_equal(res1$Std.error,res2$Std.error)
  expect_equal(res1$CI,res2$CI)
})



res4 <- cohen(m1,m2,sd1,sd2,n1,n2,conf.level,var.equal,unbiased=T, alternative)
res3 <- datacohen(Group.1,Group.2,conf.level,var.equal,unbiased=T, alternative,na.rm)

testthat::test_that("functions are consistent",{
  expect_equal(res1$Meandiff,res2$Meandiff)
  expect_equal(res1$Std.error,res2$Std.error)
  expect_equal(res1$CI,res2$CI)
})
