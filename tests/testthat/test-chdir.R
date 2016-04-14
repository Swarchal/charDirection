context("chdir")

a <- matrix(rnorm(10000), 100, 100)
a2 <- matrix(rnorm(10000), 100, 100)
aNA <- matrix(rnorm(10000), 100, 100)
aNA[1, 100] <- NA
b <- matrix(rnorm(10100), nrow = 101, ncol = 100)
a_const <- matrix(rnorm(10000), 100, 100)
a_const[1, ] <- rep(1, 100)

test_that("unequal rows returns error",{
    expect_error(chdir(a, b)) 	  
})

test_that("contains NA returns error",{
    expect_error(chdir(a, aNA))
})

test_that("constant variance returns error",{
    expect_error(chdir(a, a_const))
})
