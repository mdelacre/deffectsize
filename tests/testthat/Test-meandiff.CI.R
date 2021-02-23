res <- meandiff.CI(m1=1,m2=0,sd1=2,sd2=1.5,n1=20,n2=25, conf.level=.95, var.equal=F, alternative="two.sided")

test_that("data types correct", {
  expect_is(res,"meandiff.CI")

})

# Checking if both meandiff.CI and datameandiff.CI are consistent

Group.1 <- rnorm(10)
Group.2 <- rnorm(12)
n1 <- length(Group.1)
n2 <- length(Group.2)
m1 <- mean(Group.1)
m2 <- mean(Group.2)
sd1 <- sd(Group.1)
sd2 <- sd(Group.2)

res2 <- meandiff.CI(n1=n1,n2=n2,m1=m1,m2=m2,sd1=sd1,sd2=sd2,conf.level=.90,var.equal=T, alternative="less")
res1 <- datameandiff.CI(Group.1, Group.2, conf.level=.90, var.equal=T,alternative="less",na.rm=TRUE)

test_that("functions are consistent",{
expect_equal(res1$Meandiff,res2$Meandiff)
expect_equal(res1$Std.error,res2$Std.error)
expect_equal(res1$CI,res2$CI)
})
