# helper ----
context("colors")


# Color ----
context("color")
test_that("vec2rgb", {
    expect_equal(vec2rgb(1:3*80), "#50A0F0")
    expect_equal(vec2rgb(1:3*80/255), "#50A0F0")
    expect_equal(vec2rgb(1:3*80/255, m = 255), "#000000")
})

v <- matrix(1,3,3); v[1,] <- 0:2/3
test_that("vec2hsv", {
    expect_equal(apply(v,2,vec2hsv), diag(3)*255)
    expect_error(apply(2*v,2,vec2hsv), "values > 1 found")
})

# hue range ----
test_that("hue range", {
    expect_equal(hue.range.split(c(0,5/6), 4:2), matrix(c(0, 8,9,15,16,20)/24, nr = 2))
    #
    hrm1 <- hue.range.split(c(0, 3/6), 4:2)
    hrm2 <- hue.range.split(c(4/6, 1), 2:3)
    expect_equal(hue.range.colors(matrix(c(0,2)/6, nc = 1)), "#FFFF00")
    expect_equal(hue.range.colors(matrix(c(0,2)/6, nc = 1), only.hues = TRUE), 1/6)
    expect_equal(hue.range.colors(matrix(c(0,2)/6, nc = 1), min, only.hues = TRUE), 0)
    expect_equal(hue.range.colors(matrix(c(0,3), nc = 1)), "#FF0900")
    expect_error(hue.range.colors(matrix(c(0,300), nc = 1)), "invalid hue value")
    expect_equal(hue.range.colors(list(hrm1, hrm2)), list(c("#FF9900","#33FF00","#00FFB2"),
                                                          c("#5C00FF","#FF008A")))
})

context("rainbows")
# rainbows ----
test_that("rainbows", {
    expect_equal(hue.range.split(c(0,5/6), 4:2), matrix(c(0, 8,9,15,16,20)/24, nr = 2))
    #
    hrm1 <- hue.range.split(c(0, 3/6), 4:2)
    hrm2 <- hue.range.split(c(4/6, 1), 2:3)
    expect_equal(hue.range.colors(matrix(c(0,2)/6, nc = 1)), "#FFFF00")
    expect_equal(hue.range.colors(matrix(c(0,2)/6, nc = 1), only.hues = TRUE), 1/6)
    expect_equal(hue.range.colors(matrix(c(0,2)/6, nc = 1), min, only.hues = TRUE), 0)
    expect_equal(hue.range.colors(matrix(c(0,3), nc = 1)), "#FF0900")
    expect_error(hue.range.colors(matrix(c(0,300), nc = 1)), "invalid hue value")

    expect_equal(hue.range.colors(list(hrm1, hrm2)), list(c("#FF9900","#33FF00","#00FFB2"),
                                                          c("#5C00FF","#FF008A")))
})
