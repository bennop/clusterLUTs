#context('between') 

# all scalar
x <- 1
named.one <- 1
named.zero <- 0

names(named.one) <- x
names(named.zero) <- x

lw <- 0
hg <- 2

test_that("all scalar, fully within", {
    expect_equal(between(x, lw, hg), named.one)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.one)
    expect_equal(between(x, lw, hg, TRUE , FALSE), 1)
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(TRUE,
                                                          dimnames = list(x,
                                                                          paste(lw, hg,
                                                                                sep = '-'))))
    expect_equal(between(1, lw, hg, FALSE, FALSE), matrix(TRUE))
})

x <- 0
names(named.one) <- x

test_that("all scalar, on lower edge", {
    expect_equal(between(x, lw, hg), named.one)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.one)
    expect_equal(between(x, lw, hg, TRUE , FALSE), 1)
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(TRUE,
                                                        dimnames = list(x,
                                                                        paste(lw, hg,
                                                                              sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(TRUE))
})

x <- 2
names(named.one) <- x

test_that("all scalar, on upper edge", {
    expect_equal(between(x, lw, hg), named.one)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.one)
    expect_equal(between(x, lw, hg, TRUE , FALSE), 1)
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(TRUE,
                                                        dimnames = list(x,
                                                                        paste(lw, hg,
                                                                              sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(TRUE))
})

x <- -2
names(named.zero) <- x

test_that("all scalar, outside below", {
    expect_equal(between(x, lw, hg), named.zero)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.zero)
    expect_equal(between(x, lw, hg, TRUE , FALSE), 0)
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(FALSE,
                                                          dimnames = list(x,
                                                                          paste(lw, hg,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(FALSE))
})

x <- 3
names(named.zero) <- x

test_that("all scalar, outside above", {
    expect_equal(between(x, lw, hg), named.zero)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.zero)
    expect_equal(between(x, lw, hg, TRUE , FALSE), 0)
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(FALSE,
                                                          dimnames = list(x,
                                                                          paste(lw, hg,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(FALSE))
})

