# Color ----
context("color")
test_that("vec2rgb", {
    expect_equal(vec2rgb(1:3*80), "#50A0F0")
    expect_equal(vec2rgb(1:3*80/255), "#50A0F0")
    expect_equal(vec2rgb(1:3*80/255, m = 255), "#000000")
    expect_equal(vec2rgb(1:4*63), "#3F7EBDFC")
    expect_equal(vec2rgb(1:4*63/255), "#3F7EBDFC")
    expect_equal(vec2rgb(1:4*63/255, m = 255), "#00000000")
    expect_equal(vec2rgb(matrix(1:6*30, nrow=3)), c("#1E3C5A","#7896B4"))
    expect_equal(vec2rgb(matrix(1:8*30, nrow=4)), c("#1E3C5A78","#96B4D2F0"))
})

v <- matrix(1,3,3); v[1,] <- 0:2/3
test_that("vec2hsv", {
    expect_equal(apply(v,2,vec2hsv), diag(3)*255)
    expect_error(apply(2*v,2,vec2hsv), "values > 1 found")
})

test_that("ColorShadeRamp works",{
    expect_known_hash(ColorShadeRamp('red'), hash = "24ab5b5e67")
    expect_known_hash(ColorShadeRamp('red', space = 'rgb'), hash = "7565b58ab3")
})
test_that("color.shades",{
    expect_known_hash(color.shades(3), hash = "64e3040bfa")
    expect_known_hash(color.shades(3, c('red', 'blue')), hash = "34bad6d26d")
    expect_known_hash(color.shades(3, c('red', 'blue'), scale = 1), hash = "34bad6d26d")
    expect_known_hash(color.shades(3, c('red', 'blue'), dir = 'bright'), hash = "07cb1861c7")
    expect_known_hash(color.shades(3, c('red', 'blue'), dir = 'dark'), hash = "658831437d")
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


