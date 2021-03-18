Group.1 <- rnorm(100)
Group.2 <- rnorm(100,2)

res <- datameandiff_CI(Group.1, Group.2, conf.level=.95, var.equal=F, alternative="two.sided",na.rm=TRUE)

test_that("data types correct", {
  expect_is(res,"datameandiff_CI")

})


