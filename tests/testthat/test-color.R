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
    expect_error(vec2rgb(1:2), "wrong vector length")
    expect_error(vec2rgb(1:5), "wrong vector length")
    expect_error(vec2rgb(matrix(1, nrow = 2, ncol = 2)), "wrong matrix dimension")
    expect_error(vec2rgb(matrix(1, nrow = 5, ncol = 2)), "wrong matrix dimension")
})

test_that("vec2hsv", {
    v <- matrix(1,3,3); v[1,] <- 0:2/3
    expect_equal(apply(v,2,vec2hsv), diag(3)*255)
    expect_error(apply(2*v,2,vec2hsv), "values > 1 found")
})

test_that("color repeat", {
    rb <- matrix(0, 3, 3); rb[c(1,7)] <- 255
    expect_equal(col.rep('#000000'  , 0      ), col2rgb(NULL))
    expect_equal(col.rep('#000000'  , 1      ), matrix(0, 3, 1))
    expect_equal(col.rep('#000000'  , 2      ), matrix(0, 3, 2))
    expect_equal(col.rep('#00000080', 2      ), matrix(0, 3, 2))
    expect_equal(col.rep('#000000'  , 2, TRUE), matrix(c(0, 0, 0, 255), 4, 2))
    expect_equal(col.rep('#00000080', 2, TRUE), matrix(c(0, 0, 0, 128), 4, 2))
    expect_equal(col.rep(col2rgb('black'), 2 ), matrix(0, 3, 2))
    expect_equal(col.rep(col2rgb(c('red','black')), 1 ), rb[,1, drop = FALSE])
    expect_equal(col.rep(col2rgb(c('red','black')), 2 ), rb[,1:2])
    expect_equal(col.rep(col2rgb(c('red','black')), 3 ), rb[,1:3])
})

test_that("expand color matrix", {
    r <- matrix(0, nrow = 3, ncol = 256); r[1] <- 255; rownames(r) <- c('red','green','blue')
    ra <- col2rgb(c('red','white'), T)
    rgb.mat <- col2rgb(c('red','green','blue'))
    expect_equal(expand.colmat('red'), r)
    expect_equal(expand.colmat('red', 3), r[,1:3])
    r[4:6] <- 255
    expect_equal(expand.colmat('red', 2, 'white'), r[,1:2])
    expect_equal(expand.colmat(col2rgb('red', T), 2, 'white'), ra)
    expect_warning(expand.colmat(rgb.mat, 2), "more colors than specified output length - unchanged")
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
    hr1 <- hue.range.split(c(0,5/6), 4:2)
    r1 <- matrix(c(0, 8,9,15,16,20)/24, nrow = 2)
    hr2 <- hue.range.split(hr1, list(c(2,1,1), 3, c(1,1)))
    r2 <- matrix(c(0.00,0.1666667,0.2583333,0.375,0.6666667,0.7583333,
                   0.15,0.2416667,0.3333333,0.625,0.7416667,0.8333333),
                 nrow = 2, byrow = TRUE)

    expect_equal(hr1, r1)
    expect_equal(hr2, r2, tolerance = 1e-7)
    #
    hrm1 <- hue.range.split(c(0, 3/6), 4:2)
    hrm2 <- hue.range.split(c(4/6, 1), 2:3)
    expect_equal(hue.range.colors(matrix(c(0,2)/6, ncol = 1)), "#FFFF00")
    expect_equal(hue.range.colors(matrix(c(0,2)/6, ncol = 1), only.hues = TRUE), 1/6)
    expect_equal(hue.range.colors(matrix(c(0,2)/6, ncol = 1), min, only.hues = TRUE), 0)
    expect_equal(hue.range.colors(matrix(c(0,3)  , ncol = 1)), "#FF0900")
    expect_error(hue.range.colors(matrix(c(0,300), ncol = 1)), "invalid hue value")
    expect_equal(hue.range.colors(list(hrm1, hrm2)), list(c("#FF9900","#33FF00","#00FFB2"),
                                                          c("#5C00FF","#FF008A")))
})

test_that("hue range init", {
    m2 <- function(x) matrix(x, nrow = 2)
    # no hues -> default
    expect_equal(hue.range.init(), m2(c(0, 5/6)))
    expect_equal(hue.range.init(symm = TRUE), m2(c(0, 5/6)))
    expect_equal(hue.range.init(limits = c(.3, .8)), m2(c(.3, .8)))
    expect_equal(hue.range.init(limits = c(.3, .8), symm = TRUE), m2(c(.3, 5/6-.3)))
    ## specify hues
    hues  <- c(.2, .4, .7, 5/6)
    expect_known_hash(hue.range.init(hues), "26fc89cc09")
    expect_known_hash(hue.range.init(hues, symm = TRUE), "23ce5b42ba")
    expect_known_hash(hue.range.init(hues, limits = c(.3, .8)), "7953168a57")
    expect_warning(hue.range.init(hues, blank = 0.5), "'blank' too large")
    op <- options(verbose = TRUE)
    expect_message(hue.range.init(hues), "auto limits")
    options(op)
    #expect_known_hash(hue.range.init(hues, limits = c(.4, .8), symm = TRUE), "74bcf61d98")
    expect_known_hash(hue.range.init(hues, limits = c(.3, .8), symm = TRUE), "96106ffcfd")
    expect_warning(hue.range.init(hues, limits = c(.4, .8), symm = TRUE), "dropped hues due to limit settings: 2")

})

context("defaults")
test_that("defaults", {
    expect_known_hash(default.rgb(), "9a5206a981")
    expect_equal(default.rgb(3), c("#FF0000", "#00FF80", "#FF00FF"))
    expect_known_hash(default.hcl(), "4df6a73d55")
    expect_known_hash(default.lab(), "e7b79d9168")
})

context("rainbows")
## rainbows ----
test_that("rainbows", {
    expect_equal(hue.range.split(c(0,5/6), 4:2), matrix(c(0, 8,9,15,16,20)/24, nrow = 2))
    #
    hrm1 <- hue.range.split(c(0, 3/6), 4:2)
    hrm2 <- hue.range.split(c(4/6, 1), 2:3)
    expect_equal(hue.range.colors(matrix(c(0,2)/6, ncol = 1)), "#FFFF00")
    expect_equal(hue.range.colors(matrix(c(0,2)/6, ncol = 1), only.hues = TRUE), 1/6)
    expect_equal(hue.range.colors(matrix(c(0,2)/6, ncol = 1), min, only.hues = TRUE), 0)
    expect_equal(hue.range.colors(matrix(c(0,3)  , ncol = 1)), "#FF0900")
    expect_error(hue.range.colors(matrix(c(0,300), ncol = 1)), "invalid hue value")

    expect_equal(hue.range.colors(list(hrm1, hrm2)), list(c("#FF9900","#33FF00","#00FFB2"),
                                                          c("#5C00FF","#FF008A")))
    expect_known_hash(rainbow_lab_ramp(), "de1a6480c1")
    expect_known_hash(rainbow_lab_ramp(TRUE), "1ee4a3b0aa")
})


context("LUT file handling")
## writeLUT ----
test_that("LUT file handling", {
    tf <- tempfile()
    rl6 <- rainbow_lab(6)
    rl3 <- rl6[,1:3]
    writelut(rl6, tf)
    ##
    expect_equal(readlut(tf), matrix(rl6, nrow = 3, byrow=T))    # wrong order
    expect_equal(testthat:::safe_digest(tf), "30138132728d6414666d038b0275260c")
    #expect_known_hash(readlut(tf), hash = "d57ea18284")
    expect_equal(readlut(tf), matrix(rl6, nrow = 3, byrow = T))    # wrong order
    ##
    expect_warning(readlut(tf, length = 9), "LUT file shorter than expected: 18 bytes \\[<27\\]")
    expect_warning(readlut(tf, length = 4), "LUT file longer than expected, ignoring trailing 6 bytes")
    ##
    #expect_known_hash(readlut(tf), hash = "d57ea18284")

    ##
    expect_warning(readlut(tf, length = 9), "LUT file shorter than expected: 18 bytes \\[<27\\]")
    expect_warning(readlut(tf, length = 4), "LUT file longer than expected, ignoring trailing 6 bytes")
    ##
    expect_equal(colmat2lutfile(rl6, tf), rl6)
    expect_equal(testthat:::safe_digest(tf), "56db3b781f0c6c5703bc2b9925205b5c")
    expect_equal(readlut(tf), rl6)                        # right order
    expect_known_hash(colmat2lutfile(rl3, tf, fill = TRUE   ), "4e9f7c64c8")
    expect_known_hash(colmat2lutfile(rl3, tf, fill = TRUE, 2), "9b5c898cf1")
    expect_known_hash(colmat2lutfile(rl3, tf, fill = TRUE, 3), "9b5c898cf1") # same!
    expect_known_hash(colmat2lutfile(rl3, tf, fill = TRUE, 4), "6f3116d02b")
    expect_known_hash(colmat2lutfile(rl3, tf, fill = TRUE, 5), "1435cc0c9c")
    expect_known_hash(colmat2lutfile(rl3, tf, fill = TRUE, 8), "4063b58304")
    expect_known_hash(colmat2lutfile(rl3, tf, fill = TRUE, 2, bw=FALSE), "9b5c898cf1")
    expect_known_hash(colmat2lutfile(rl3, tf, fill = TRUE, 3, bw=FALSE), "9b5c898cf1")
    expect_known_hash(colmat2lutfile(rl3, tf, fill = TRUE, 4, bw=FALSE), "c3d1e4e3a7")
    expect_known_hash(colmat2lutfile(rl3, tf, fill = TRUE, 8, bw=FALSE), "348d49c33a")


    ##
    rl6[1,1] <- 256
    expect_equal(writelut(rl6, tf), 99)
    expect_equal(file.exists(tf), FALSE)   # file should have been deleted
    expect_error(colmat2lutfile(rl6, tf), "error writing LUT")
    ##
    ## not multiple of 3 bytes (8 bytes)
    writeBin(pi, tf)
    suppressWarnings( expect_known_hash(readlut(tf), "7878c28d1e") )
    expect_error(readlut(tf, 3), "^LUT.*LUT$")
 })
