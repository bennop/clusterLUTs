context("color")
## Color ----
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

test_that("vec2hsv", {
    v <- matrix(1,3,3); v[1,] <- 0:2/3
    expect_equal(apply(v,2,vec2hsv), diag(3)*255)
    expect_error(apply(2*v,2,vec2hsv), "values > 1 found")
})

test_that("ColorShadeRamp works",{
    expect_known_hash(ColorShadeRamp('red'), hash = "24ab5b5e67")
    expect_known_hash(ColorShadeRamp('red', space = 'rgb'), hash = "7565b58ab3")
})

test_that("color.shades",{
    c.rb <- c('red', 'blue')
    c.rgb <- c('red', 'green', 'blue')
    expect_known_hash(color.shades(1), hash = "fe96ed264e")
    expect_known_hash(color.shades(3), hash = "64e3040bfa")
    expect_known_hash(color.shades(3, c.rb), hash = "34bad6d26d")
    expect_known_hash(color.shades(2:3, c.rgb), hash = "7cf71e94d7")
    expect_known_hash(color.shades(3, c.rb, scale = 1), hash = "34bad6d26d")
    expect_known_hash(color.shades(3, c.rb, dir = 'bright'), hash = "07cb1861c7")
    expect_known_hash(color.shades(3, c.rb, dir = 'dark'), hash = "658831437d")
    expect_known_hash(color.shades(3:4), hash = "8cb5c08b19")
    expect_known_hash(color.shades(3:4, c.rb), hash = "9245e933f9")
    # gscale = FALSE
    expect_known_hash(color.shades(1:3, c.rgb, gscale = FALSE), hash = "f48b4f796a")
    expect_known_hash(color.shades(2:4, c.rgb, dir = 'd', gscale = FALSE), hash = "794b3f5fc6")
    expect_known_hash(color.shades(2:4, c.rgb, dir = 'b', gscale = FALSE), hash = "b545e28063")
    expect_known_hash(color.shades(3:4, gscale = FALSE), hash = "8f658d35aa")
    expect_known_hash(color.shades(3:4, c.rb, gscale = FALSE), hash = "7a0e7708d0")
})

test_that("cutshades",{
    set.seed(42)
    expect_known_hash(cutshades(ct <-  sample(6, 15, TRUE)), hash = "96125ec894")
    expect_known_hash(color.shades(3, c('red', 'blue'), scale = 1), hash = "34bad6d26d")
    expect_known_hash(color.shades(3, c('red', 'blue'), dir = 'bright'), hash = "07cb1861c7")
    expect_known_hash(color.shades(3, c('red', 'blue'), dir = 'dark'), hash = "658831437d")
})


## hue range ----
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
## rainbows ----
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


context("LUT file handling")
## writeLUT ----
test_that("LUT file handling", {
    tf <- tempfile()
    rl6 <- rainbow_lab(6)
    writelut(rl6, tf)
    ##
    expect_equal(testthat:::safe_digest(tf), "30138132728d6414666d038b0275260c")
    expect_known_hash(readlut(tf), hash = "d57ea18284")
    ##
    expect_warning(readlut(tf, length = 9), "LUT file shorter than expected: 18 bytes \\[<27\\]")
    expect_warning(readlut(tf, length = 4), "LUT file longer than expected, ignoring trailing 6 bytes")
    ##
    rl6[1,1] <- 256
    expect_equal(writelut(rl6, tf), 99)
    expect_equal(file.exists(tf), FALSE)   # file should have been deleted
    ##
    ## not multiple of 3 bytes
    writeBin(pi, tf)
    suppressWarnings( expect_known_hash(readlut(tf), "7878c28d1e") )
    expect_error(readlut(tf, 3), "^LUT.*LUT$")
 })
