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
    t1 <- table(sample(1:2, 10, T))
    t2 <- table(sample(7:9, 10, T))
    m1 <- matrix(c(1:2, 5,5), nrow=2, byrow = TRUE)
    m2 <- matrix(c(7:9, 3,3,4), nrow=2, byrow = TRUE)
    m2n <- m2; rownames(m2n) <- c('downlink', 'n')
    tl1 <- list(t1, t2)
    tl2 <- list(t1, m2)
    tl2b <- list(m2, t1)

    tl3 <- list(t1, m2n)
    tl4 <- list(m1, m2)
    tl5 <- list(m1, t(m2))
    res0  <- matrix(c(1,2,7,8,9, 5,5,3,3,4), nr = 2, byrow = T,
                    dimnames = list(c('downlink', 'n'), NULL))
    res1  <- matrix(c(1,1,2,2,2, 1,2,7,8,9, 5,5,3,3,4), nr = 3, byrow = T,
                   dimnames = list(c('index', 'downlink', 'n'), NULL))
    res1b <- matrix(c(1,1,1,2,2, 7,8,9,1,2, 3,3,4,5,5), nr = 3, byrow = T,
                    dimnames = list(c('index', 'downlink', 'n'), NULL))
    expect_equal(concat.tbl.list(tl1), res0)
    expect_equal(concat.tbl.list(tl1, idx = TRUE), res1)
    expect_equal(concat.tbl.list(tl2, idx = TRUE), res1)
    expect_equal(concat.tbl.list(tl2b, idx = TRUE), res1b)
    expect_equal(concat.tbl.list(tl3, idx = TRUE), res1)
    expect_equal(concat.tbl.list(tl4, idx = TRUE), res1)
    expect_error(concat.tbl.list(tl5, idx = TRUE), "incompatible inputs")

})

test_that("table 2 matrix",{
    set.seed(42)
    tbl  <- table(sample(1:3, 10, TRUE))
    exm1 <- matrix(0, 2,3)
    exm1b <- exm1; rownames(exm1b) <- letters[1:2]
    exm2 <- matrix(0, 3,2)
    rm1 <- matrix(c(1,5,2,3,3,2), nrow=2, dimnames = list(c('idx','n'), 1:3))
    rm2 <- rbind(rm1, 4:6)
    rm2b <- rm2; rownames(rm2b)[3] <- 'uplink'
    rm3 <- rbind(rm1, exm1)

    expect_equal(tbl2mat(tbl), rm1)
    expect_equal(tbl2mat(tbl, 4:6), rm2b)
    expect_equal(tbl2mat(tbl, 4:6, vname = 'ext'), rm2b)
    expect_warning(tbl2mat(tbl, 4:5), "incompatible vector 'uplink' - ignored")
    expect_equal(tbl2mat(tbl, exm1), rm3)
    expect_warning(tbl2mat(tbl, matrix(0,2,2)), "incompatible matrix 'uplink' - ignored")
    expect_warning(tbl2mat(tbl, array(1, dim = rep(1,3))), "incompatible array 'uplink' - ignored")

})

test_that("re-index",{
    ct <- c(4, 1, 2, 3, 1, 1, 1, 1, 4, 1, 1, 3, 1, 3)
    res  <- c(13, 1, 9,10, 2, 3, 4, 5,14, 6, 7,11, 8,12)
    names(res) <- c(4.1,1.1,2.1,3.1,1.2,1.3,1.4,1.5,4.2,1.6,1.7,3.2,1.8,3.3)
    expect_equal(reidx.cut(ct), res)
})

# tree ----
context("tree")
test_that("read tree", {
    expect_known_hash(read.tree(), hash = "23ef87605d")
    expect_known_hash(read.tree('sbm_1_145_0_atlas.mat',
                                system.file('extdata',
                                            package = 'clusterLUTs')), hash = "23ef87605d")
    tf <- tempfile()
    expect_error(read.tree(tf), "^tree file not found")
})

test_that("tree.ranges", {
    expect_known_hash(tree.ranges(dummy.tree()), "80d848e629")
})

test_that('subtables', {
    expect_known_hash(subtables(tree.ranges(dummy.tree())), "16945ffa45")
    expect_error(subtables(letters), "no sub.tables found in ")

})

test_that("treeluts", {
    cat("   preliminary\n")
    expect_known_hash(treeluts(dummy.tree()), "e89e26de0c")
})


# %&% ----
context("%&%")
test_that("&", {
    expect_equal(7 %&% 5, 5)
})

# debug ----
context("debug")
test_that("dummy tree", {
    expect_known_hash(dummy.tree(  ), hash = "d4032aea3f")
    expect_known_hash(dummy.tree(25), hash = "f863368944")
})

test_that("randomize dt", {
    set.seed(42)
    expect_known_hash(randomize.cutree(dummy.tree()), "469fc7b26c")
})

test_that("dend.with.cuts", {
    dt <- dummy.tree()
    expect_equal(dend.with.cuts(dt),
                 c(44.89947, 43.97886, 42.27916), tolerance = .00002)
    expect_known_hash(dend.with.cuts(dt), "1b04029414")
    expect_known_hash(dend.with.cuts(dt, cuts = 4:6), "5c77f90029")
    expect_error(dend.with.cuts(matrix(0,1,1)),
                 "not recognized as either hclust or dummy.tree output")
})

test_that("hri.plot", {
    expect_known_hash(clusterLUTs:::hri.plot(1:5/6), "d5a90bb13d")
})
