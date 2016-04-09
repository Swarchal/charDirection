context("shrink_mat")

# example data
data(iris)
R <- nipals(iris[,1:4], 3)$T

test_that("returns a matrix",{
    tmp <- shrink_mat(R, samples_count = 10, r = 1)
    expect_is(tmp, "matrix")
    expect_true(is.numeric(tmp))
})

test_that("samples_count affects matrix",{
    tmp1 <- shrink_mat(R, samples_count = 3, r = 1)
    tmp2 <- shrink_mat(R, samples_count = 100, r = 1)
    expect_true(any(tmp1 != tmp2))
})

test_that("r affects matrix",{
    tmp1 <- shrink_mat(R, samples_count = 50, r = 1)
    tmp0 <- shrink_mat(R, samples_count = 50, r = 0)
    var1 <- var(diag(tmp1))
    var2 <- var(diag(tmp0))
    expect_true(var1 > var2)
})
