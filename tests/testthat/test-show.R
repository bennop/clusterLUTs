test_that("show brain LUT", {
    expect_equal(show.brain.lut(1, matrix(0,3,4), matrix(c(1,2,1,1),2), test = TRUE), "direct lut")
    expect_equal(show.brain.lut(1, matrix(0,2,3), matrix(c(1,2,1,1),2), test = TRUE), "t(lut)")
    #    expect_failure(show.brain.lut(1, matrix(0,2,2), matrix(c(1,2,1,1),2)))
})
