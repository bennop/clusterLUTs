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
dt <- dummy.tree()
test_that("dummy tree", {
    expect_known_hash(dt, hash = "d4032aea3f")
    expect_equal(dend.with.cuts(attr(dt, 'hc'), as.numeric(colnames(dt)), col = "#ffa050"),
                  c(44.89947, 43.97886, 42.27916), tolerance = .00002)
})
