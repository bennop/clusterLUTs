# helper ----
context("helper")

m <- matrix(1:6, ncol = 2)
test_that("matrix to list", {
    expect_equal(mat2list(m), list(1:3, 4:6))
    expect_equal(mat2list(m, 1), list(c(1,4), c(2,5), c(3,6)))
})

test_that("vlevels", {
    expect_equal(vlevels(c(1:4, 1:3, 1:2, 1)), 4)
})

test_that("concat.tbl.list", {
    set.seed(42)
    tl <- list(table(sample(1:2, 10, T)), table(sample(7:9, 10, T)))
    res  <- matrix(c(1,1,2,2,2, 1,2,7,8,9, 5,5,3,3,4), nr = 3, byrow = T,
                   dimnames = list(c('index', 'downlink', 'n'), NULL))
    expect_equal(concat.tbl.list(tl), res[2:3,])
    expect_equal(concat.tbl.list(tl, idx = TRUE), res)
})

# color ----
context("color")
test_that("vec2rgb", {
    expect_equal(vec2rgb(1:3*80), "#50A0F0")
    expect_equal(vec2rgb(1:3*80/255), "#50A0F0")
    expect_equal(vec2rgb(1:3*80/255, m = 255), "#000000")
})

v <- matrix(1,3,3); v[1,] <- 0:2/3
test_that("vec2hsv", {
    expect_equal(apply(v,2,vec2hsv), diag(3)*255)
})

# hue range ----
test_that("hue range", {
    expect_equal(split.hue.range(c(0,5/6), 4:2), matrix(c(0, 8,9,15,16,20)/24, nr = 2))
    expect_equal(hue.range.colors(matrix(c(0,2)/6, nc = 1)), "#FFFF00")
    expect_equal(hue.range.colors(matrix(c(0,2)/6, nc = 1), only.hues = TRUE), 1/6)
    expect_equal(hue.range.colors(matrix(c(0,2)/6, nc = 1), min, only.hues = TRUE), 0)
    expect_error(hue.range.colors(matrix(c(0,3), nc = 1)), "invalid hsv color")
})

# debug ----
context("debug")
test_that("dummy tree", {
    dt <- dummy.tree()
    expect_known_hash(dt, hash = "d4032aea3f")
    expect_equal(dend.with.cuts(dt, cut.col = "#ffa050"),
                  c(44.89947, 43.97886, 42.27916), tolerance = .00002)
})
