#context('between')
#
# all scalar  ----

context("between - all scalar")
x <- 1
named.one <- 1
named.zero <- 0

names(named.one) <- x
names(named.zero) <- x

lw <-  1
test_that("between - scalar, pointwise",{
    expect_equal(between(x, lw           ), named.one)
    # skippng high while not naming the other parameters fails
    expect_error(between(x, lw, T, T, 'b'), "apparently no `high' but unnamed parameters that take it's place - please specify", fixed = TRUE)
    # either specify ALL later parameters that could slide in ...
        expect_equal(between(x, lw, i=T, n=T, e='b'), named.one)
    # ... or fill in high
    expect_equal(between(x, lw, NULL, T, T, 'b'), named.one)
    expect_equal(between(x, lw, i=T, n=T, e='n'), named.one)
    expect_equal(between(x, lw, i=T, n=T, e='l'), named.one)
    expect_equal(between(x, lw, i=T, n=T, e='r'), named.one)
    expect_equal(between(x, lw, i=T, n=F, e='b'), 1)
    expect_equal(between(x, lw, i=T, n=F, e='n'), 1)
    expect_equal(between(x, lw, i=T, n=F, e='l'), 1)
    expect_equal(between(x, lw, i=T, n=F, e='r'), 1)
    expect_equal(between(x, lw, i=F, n=T, e='b'), matrix(TRUE, dimnames = list(x, paste(lw, lw,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, i=F, n=T, e='n'), matrix(TRUE, dimnames = list(x, paste(lw, lw,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, i=F, n=T, e='l'), matrix(TRUE, dimnames = list(x, paste(lw, lw,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, i=F, n=T, e='r'), matrix(TRUE, dimnames = list(x, paste(lw, lw,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, i=F, n=F, e='b'), matrix(TRUE))
    expect_equal(between(x, lw, i=F, n=F, e='n'), matrix(TRUE))
    expect_equal(between(x, lw, i=F, n=F, e='l'), matrix(TRUE))
    expect_equal(between(x, lw, i=F, n=F, e='r'), matrix(TRUE))
})

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
names(named.zero) <- x

test_that("all scalar, on lower edge", {
    expect_equal(between(x, lw, hg), named.one)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.one)
    expect_equal(between(x, lw, hg, TRUE , FALSE), 1)
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(TRUE,
                                                          dimnames = list(x,
                                                                          paste(lw, hg,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(TRUE))
    expect_equal(between(x, lw, hg, e = 'r'), named.zero)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE , e = 'r'), named.zero)
    expect_equal(between(x, lw, hg, TRUE , FALSE, e = 'r'), 0)
    expect_equal(between(x, lw, hg, FALSE, TRUE , e = 'r'), matrix(FALSE,
                                                                   dimnames = list(x,
                                                                                   paste(lw, hg,
                                                                                         sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE, e = 'r'), matrix(FALSE))
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

# vector / scalar ----

context("between - vector x, scalar range")
x <- 1:2 - 0.5
lx <- length(x)
named.one <- rep(1, length(x))
named.zero <- rep(0, length(x))

names(named.one) <- x
names(named.zero) <- x

lw <- 0
hg <- 2
ll <- length(lw)

test_that("vector x, scalar range, fully within", {
    expect_equal(between(x, lw, hg), named.one)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.one)
    expect_equal(between(x, lw, hg, TRUE , FALSE), rep(1, lx))
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(rep(TRUE, lx*ll),
                                                          dimnames = list(x,
                                                                          paste(lw, hg,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(rep(TRUE, lx*ll), nc = ll))
})

test_that("vector x, scalar range, fully within", {
    x <- 0:1
    names(named.one) <- x

    expect_equal(between(x, lw, hg), named.one)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.one)
    expect_equal(between(x, lw, hg, TRUE , FALSE), rep(1, lx))
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(rep(TRUE, lx*ll),
                                                          dimnames = list(x,
                                                                          paste(lw, hg,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(rep(TRUE, lx*ll), nc = ll))
})

test_that("vector x, scalar range, inside / on upper edge", {
    x <- 1:2
    names(named.one) <- x

    expect_equal(between(x, lw, hg), named.one)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.one)
    expect_equal(between(x, lw, hg, TRUE , FALSE), rep(1, lx))
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(rep(TRUE, lx*ll),
                                                          dimnames = list(x,
                                                                          paste(lw, hg,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(rep(TRUE, lx*ll), nc = ll))
})


test_that("vector x, scalar range, on upper edge / outside", {
    x <- 2:3
    names(named.one) <- x
    oz = 1:0; names(oz) <- x
    expect_equal(between(x, lw, hg), oz)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), oz)
    expect_equal(between(x, lw, hg, TRUE , FALSE), 1:0)
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(c(TRUE,FALSE), nc=1,
                                                          dimnames = list(x,
                                                                          paste(lw, hg,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(c(TRUE,FALSE), nc=1))
})

test_that("vector x, scalar range, outside both sides", {
    x <- c(-1,3)
    named.zero <- rep(0,lx)
    names(named.zero) <- x

    expect_equal(between(x, lw, hg), named.zero)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.zero)
    expect_equal(between(x, lw, hg, TRUE , FALSE), c(0,0))
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(rep(FALSE,2), nc = 1,
                                                          dimnames = list(x,
                                                                          paste(lw, hg,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(rep(FALSE,2), nc = 1))
})

x <- 3:4
names(named.zero) <- x

test_that("vector x, scalar range, outside above", {
    expect_equal(between(x, lw, hg), named.zero)    # default, same as next
    expect_equal(between(x, lw, hg, TRUE , TRUE ), named.zero)
    expect_equal(between(x, lw, hg, TRUE , FALSE), c(0,0))
    expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(rep(FALSE,2), nc = 1,
                                                          dimnames = list(x,
                                                                          paste(lw, hg,
                                                                                sep = '-'))))
    expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(rep(FALSE,2), nc = 1))
})

# scalar / vector ----
context("between - scalar x / vector range")
lw <- c(0,2)
hg <- c(1,3)


test_that("scalar x, vector range, outside above", {
    # expect_equal(between(x, lw, hg), named.zero)    # default, same as next
    # expect_equal(between(x, lw, hg, TRUE , TRUE ), named.zero)
    expect_equal(between(-1, lw, hg, TRUE , FALSE), 0)
    # expect_equal(between(x, lw, hg, FALSE, TRUE ), matrix(rep(FALSE,2), nc = 1,
    #                                                       dimnames = list(x,
    #                                                                       paste(lw, hg,
    #                                                                             sep = '-'))))
    # expect_equal(between(x, lw, hg, FALSE, FALSE), matrix(rep(FALSE,2), nc = 1))
})

# all vector
