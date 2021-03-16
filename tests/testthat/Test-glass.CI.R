res <- glass.CI(m1,m2,sd1,sd2,n1,n2,conf.level,unbiased=T,alternative="two.sided")

test_that("data types correct", {
  expect_is(res,"glass.CI")

})

# Checking if both glass.CI and dataglass.CI are consistent

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

res2 <- glass.CI(m1,m2,sd1,sd2,n1,n2,conf.level,unbiased=F, alternative)
res1 <- dataglass.CI(Group.1,Group.2,conf.level,unbiased=F, alternative,na.rm)

testthat::test_that("functions are consistent",{
  expect_equal(res1$Meandiff,res2$Meandiff)
  expect_equal(res1$Std.error,res2$Std.error)
  expect_equal(res1$CI,res2$CI)
})



res4 <- glass.CI(m1,m2,sd1,sd2,n1,n2,conf.level,unbiased=T, alternative)
res3 <- dataglass.CI(Group.1,Group.2,conf.level,unbiased=T, alternative,na.rm)

testthat::test_that("functions are consistent",{
  expect_equal(res4$Meandiff,res3$Meandiff)
  expect_equal(res4$Std.error,res3$Std.error)
  expect_equal(res4$CI,res3$CI)
})



