context("chdir")

a <- matrix(rnorm(10000), 100, 100)
a2 <- matrix(rnorm(10000), 100, 100)
aNA <- matrix(rnorm(10000), 100, 100)
aNA[1, 100] <- NA
b <- matrix(rnorm(10100), nrow = 101, ncol = 100)
a_const <- matrix(rnorm(10000), 100, 100)
a_const[1, ] <- rep(1, 100)

test_that("unequal rows returns error",{
    expect_error(chdir(a, b, samples = 1:100)) 	  
})

test_that("contains NA returns error",{
    expect_error(chdir(a, aNA, samples = 1:100))
})

test_that("constant variance returns error",{
    expect_error(chdir(a, a_const, samples = 1:100))
    expect_error(chdir(a_const, a, samples = 1:100))
})

test_that("runs without error",{
    expect_is(chdir(a, a2, samples = 1:100), 'matrix')
})
