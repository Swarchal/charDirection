context("nipals")

# example data
data(iris)

test_that("returns a list",{
    out <- nipals(iris[, 1:4], 2)
    expect_is(out, "list")
    expect_equal(length(out), 4)
    expect_equal(names(out), c("T", "P", "pcvar", "precision"))
})

test_that("returns expected number of principal components",{
    out1 <- nipals(iris[, 1:4], 1)
    out3 <- nipals(iris[, 1:4], 3)
    expect_equal(ncol(out1$T), 1)
    expect_equal(ncol(out3$T), 3)
})

test_that("errors when more components than initial columns",{
    expect_error(nipals(iris[, 1:4], 5))
})

