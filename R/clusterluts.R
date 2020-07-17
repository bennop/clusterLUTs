
#' local project directory
#'
#' Used as default for \code{\link{tree.luts}} and \code{\link{read.tree}},
#' simply check for existence and return input (or default).
#'
#'
#' ~~~~~~~~~~~~~~~~~~~ C A V E A T ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#'       Adjust local defaults
#'
#'  these settings define DEFAULTS that depend on my directory structure and
#'  need to be adjusted to work for others
#'
#' @param basedir project directory
#' @param ... ignored
#'
#' @return path to project dirctory
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' projdir()
projdir <- function(basedir = '~/Work/4philipp/BrainLUTs',
                    ...){
    dir.not.found <- function(dir){
        sprintf("%d dir not found - please specify explicitely", dir)
    }
    if(!dir.exists(basedir)) stop(dir.not.found('base'))

     return(basedir)

}

#' local LUT directory
#'
#' Used as default for \code{\link{tree.luts}}
#'
#'
#' ~~~~~~~~~~~~~~~~~~~ C A V E A T ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#'       Adjust local defaults
#'
#'  these settings define DEFAULTS that depend on my directory structure and
#'  need to be adjusted to work for others
#'
#' @param basedir parent directory, usually project
#' @param subdir LUT directory
#' @param ... ignored
#'
#' @return path to LUT dirctory
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' LUTdir()
LUTdir <- function(basedir = projdir(),
                   subdir = 'luts',
                   ...){
    dir.not.found <- function(dir){
        sprintf("%d dir not found - please specify explicitely", dir)
    }
    if(!dir.exists(basedir)) stop(dir.not.found('base'))

    LUTdir <- file.path(basedir, subdir)
    if(!dir.exists(LUTdir)) stop(dir.not.found('LUT'))
    return(LUTdir)

}


## needed packages ----
require(R.matlab)    # for readMat()
require(gtools)      # for even()

##---- Rcpp helper ----------------
## source / load Rcpp helper function
## fix path problem later
# Rcpp::sourceCpp("writelut.cpp")

##---  R functions -----------------
#' LUT generation for hierachical clustering tree
#'
#' Main function of this package.
#' \code{clusters} defines the cluster assignments per cut. The number of clusters
#' is taken from the unique values per column (need not be regular).
#'
#' In order to show the substructure of clusters ... (\code{\link{cut.shades}})
#'
#' The LUTs are named \code{sprintf("%s%03d.lut",basename,n)} where \code{n} is
#' the number of clusters in the cut
#'
#' @param clusters matrix with cluster assignments, one column per cut
#' @param outdir where to write the LUTs to
#' @param basename prefix for LUT names ['lut']
#' @param ... ignored
#'
#' @return
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' treeluts(cbind(sample(1:10,20,TRUE), 1:20), '.', "example")
treeluts <- function(clusters = read.tree(),
                     outdir   = LUTdir(),
                     basename = 'lut',
                     ...){
    nc <- nrow(clusters)
    lut.length <- 256     # required by MRIcron
    last.cut <- ncol(clusters)
    apply(clusters, 2,   # for each cut
          function(v){
              n <- vlevels(v)

              #pal <- c(rainbow(n, end = 0.65))
              #fullpal <- pal[v]
              fullpal <- cut.shades(v, ...)
              pmat <- c(fullpal,                      # colors
                        rep(col2rgb('#000000'),       # more black
                            lut.length - nc - 1),
                        col2rgb(ifelse(nc<lut.length, # extra slot?
                                       '#FFFFFF',     # end with white
                                       NULL)
                        )          # don't
              )

              imat <- as.integer(as.vector(t(pmat)))
              imat[imat>128] <- imat[imat>128] - 256
              file <- path.expand(file.path(outdir,
                                            sprintf("%s%03d.lut",
                                                    basename,
                                                    n)))
              cat(file,"\n")
              writelut(imat,
                       file)
              return(n)
          }
    )
}

#' read lookup table
#'
#' The file will be read as binary with three (unsigned) bytes per color and
#' return a color matrix with columns corresponding to the color entries (RGB).
#'
#' If \code{length} is specified, read only that many entries, otherwise all.
#' Some sonsistency checks on the file size are performed.
#'
#' @param file lookup table file
#' @param length number of color entries to read
#' @param ... ignored
#'
#' @return color matrix (see description)
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' \dontrun{
#'    lut <-  readlut('test.lut')
#' }
readlut <- function(file,
                    length = file.size(file)/3,
                    ...){
    fsize <- file.size(file)
    if(fsize < 3 * length){
        if (fsize %% 3 == 0){
            warning("LUT file shorter than expected: ", fsize, " bytes [<", 3*length, "]")
            length <- fsize/3
        } else {
            stop("LUT file size (", fsize, ") not multiple of 3, thus likely not a LUT")
        }
    }
    l <- readBin(file,
                 "integer",
                 size   = 1,
                 signed = FALSE,
                 n      = 3 * length)
    return(matrix(l,
                  nrow  = 3,
                  byrow = F))
}

#' read cluster tree information as obtained through Matlab
#'
#' The mat-file is expected to have one component "cluster.info" which, in turn
#' contains a list of lists. Each of the latter is a \eqn{nx1} matrix with cluster
#' assignments (the number of cluters at each cut is given by the number of unique
#' entries)
#'
#' Each of these matrices (vectors) is forms one column in the output matrix as
#' would be obtained by \code{\link[stats]{cutree}}. The matrix is given column
#' names indicating the number of clusters for that cut.
#'
#' @param file Matlab hclust result file (Matlab *.mat-file)
#'
#' @return a matrix with just the cluster assignments
#' @export
#' @importFrom R.matlab readMat
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' \dontrun{
#'    tree <- read.tree(path_to_treefile)
#' }
read.tree <- function(file = file.path(projdir(),
                                       'bm_1_145_0_result_hclust_atlas.mat')){
    if(!file.exists(file)){
        stop(sprintf("tree file not found\n\t'%s'",
                     file))
    }
    ci <- readMat(file)$cluster.info
    tree <- sapply(ci, function(l)l[[1]])
    cuts <- apply(tree, 2, vlevels)
    colnames(tree) <- cuts
    return(invisible(tree))
}

#' Create a color ramp with  shades
#'
#' The shades will range from black over \code{col} to white.
#' for symmetry the color should be a color with full saturation
#' @param col base color
#' @param space 'Lab' or 'rgb', passed on to \code{\link[grDevices]{colorRamp}}
#' @param ... also  passed on to \code{\link[grDevices]{colorRamp}}
#'
#' @return a function mapping [0,1] to (black ... \code{col} ... white)
#'
#' @importFrom grDevices colorRamp
#' @export
#'
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' ColorShadeRamp('red')
ColorShadeRamp <- function(col,
                           space = 'Lab',
                           ...){
    colorRamp(c('black',
                col,
                'white'),
              space = space,
              ...)
}

#' create color shades
#'
#' For each color in \code{col} create as many shades as indicated by the value
#' of \code{reps} at the corresponding index by grading towards 'white' - the
#' shades are going \code{scale} of the way to white.
#'
#' If \code{reps} is scalar it is recycled, otherwise it is filled with ones should
#' it be shorter than \code{col}.
#'
#' The shades can be constructed either towards white, black, or, symmetrically.
#' This is determined by setting \code{direction} to 'bright', 'dark', or 'all',
#' respectively. The latter is the default.
#'
#' @param col color vector; must be a valid argument to \code{\link[grDevices]{col2rgb}}
#'  or a \eqn{3xN} RGB color matrix (one column per color).
#' @param reps vector of repetitions
#' @param scale how much of the range to white should be covered
#' @param direction which way to build the shades, see Description
#' @param ... passed to \code{\link{colorShadeRamp}}
#'
#' @return Color matrix with on entry per column
#' @importFrom grDevices col2rgb
#' @importFrom gtools even
#' @export
#'
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' color.shades(c('red', 'blue'), 3)
#' show.shades(color.shades(c('red', 'green', 'blue'), 2:4))
color.shades <- function(col,
                         reps,
                         scale = ifelse(direction == 'bright', 0.7, 0.5),
                         direction = c('all', 'bright', 'dark'),
                         ...){
    # handle parameters
    direction <- match.arg(direction)
    scale <- min(0.9, max(0.1, scale))

    #
    if(is.null(dim(col))){
        col <- col2rgb(col)
    }
    len.diff <- ncol(col) - length(reps)
    if(len.diff > 0){
        reps <- if(length(reps == 1)){   # scalar:
            rep(reps, ncol(col))       #    replicate
        } else {                         # vector:
            c(reps, rep(1, len.diff))    #    use 1 for positions not specified
        }
    }

    col.reps <- list()

    for (cidx in 1:ncol(col)){
        # color function (black - color - white)
        csf <- ColorShadeRamp(vec2rgb(col[,cidx]), ...)

        # the max() contruct is needed to avoid division by
        # zero when reps[cidx]==1
        rmp.n <- ifelse(direction == 'all' && even(reps[cidx]),
                        reps[cidx]+1 ,
                        reps[cidx])

        ramp.vals <- if(rmp.n == 1){
            0.5
        } else {
            switch(direction,
                   all = {
                       v <- seq(0.25,
                                0.75,
                                length.out = rmp.n)
                       if(even(reps[cidx])){
                           v[-1]
                       } else {
                           v
                       }
                   },
                   bright = seq(0.5, 1-(1-scale)*0.5,
                                length.out = rmp.n),
                   dark =  seq(0.5, (1-scale)*0.5,
                               length.out = rmp.n)
            )
        }
        col.reps[[cidx]] <-  t(csf(ramp.vals))
    }
    shades <- matrix(unlist(col.reps), nrow = 3)
    attr(shades, 'reps') <- reps
    attr(shades, 'dir') <- direction
    return(shades)
}



#' re-index cut(ree) vector
#'
#' internal rank increasing within ties so that after reordering the order of ties
#' can be preseerved
#'
#' the output vector is named with 1000 * original level + position with ties for
#' that level, see example
#'
#' This is useful to reorder, e.g., a color LUT from \code{\link{color.shades}}
#' according to levels.
#'
#' @param cut cut from hierarchical clustering
#'
#' @return rank vector with names to indicate level and order within ties
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' cut <- c(4, 1, 2, 3, 1, 1, 1, 1, 4, 1, 1, 3, 1, 3)
#' reidx.cut(cut)
#' ## 4001 1001 2001 3001 1002 1003 1004 1005 4002 1006 1007 3002 1008 3003
#' ##   13    1    9   10    2    3    4    5   14    6    7   11    8   12
#' ## the first occurence of level 1 is named 1001 (see details)
#' cs <- color.shades(rainbow(4), table(t1))
#' rocs <- cs[, reidx.cut(t1)]
reidx.cut <- function(cut){
    tbl <- table(cut)
    vals <- 1:length(unique(cut))
    for (i in vals)
        cut[cut == i] <- 1000*i + 1:tbl[i]
    rcut <- rank(cut)
    names(rcut) <- cut
    return(rcut)
}



#' RGB-based rainbow
#'
#' @param n number of colors
#' @param ...  passed to rainbow_rgb
#'
#' @return rainbow colors
#' @importFrom grDevices rainbow
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' default.rgb()
default.rgb <- function(n = 12, ...){
    rainbow(n, end = 5/6)
}

#' HCL-based rainbow
#'
#' @param n number of colors
#' @param ... passed to rainbow_hcl
#'
#' @return rainbow colors
#' @importFrom colorspace rainbow_hcl
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' default.hcl()
default.hcl <- function(n = 12, ...){
    rainbow_hcl(n, end = 5/6)
}

#' Lab-based rainbow
#'
#' @param n number of colors
#' @param ... passed to rainbow_lab
#'
#' @return rainbow colors
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' default.lab()
default.lab <- function(n = 12, ...){
    rainbow_lab(n, end = 5/6, ...)
}

#' prepare shades and reorder to match cut
#'
#' For the given cut \code{cut} analyze the internal structure of each cluster
#' and use that as basis for a shading through \code{\link{color.shades}}. Then
#' reorder the shading to arrange according to cluster members.
#'
#'
#'
#' @param cut vector of cluster assignments
#' @param col.fun function mapping a number to that many colors from a palette
#'
#' @return color shades reordered (see details)
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' cut.shades(sample(6, 15, T))
cut.shades <- function(cut,
                       col.fun = default.rgb){
    tbl <- table(cut)
    cs <- color.shades(col.fun(length(tbl)),
                       tbl)
    cso <- cs[, order(reidx.cut(cut))]
    attr(cso, 'reps') <- attr(cs, 'reps')
    attr(cso, 'cut') <- cut
    return(cso)
}

##----  Helper functions, mostly for internal use ----

#' between test
#'
#' test for \code{x} to lie between \code{low} and \code{high} (including the
#' margins), i.e., \eqn{low \leq x \leq high}
#'
#' \code{between} is vectorized for \code{x} as well as for the limits. In the
#' case of vectorized limits \code{low} and \code{high} have to have the same
#' length and the ranges are defined by corresponding pairs from \code{low} and
#' \code{high} (those ranges are allowed to overlap).
#'
#'
#' @param x     value(s) to test
#' @param low   lower bound(s) of test interval(s)
#' @param high  upper bound(s) of test interval(s)
#' @param index when set return indices, otherwise match matrix (see description)
#' @param named whether to put names on the return value
#' @param ... ignored
#'
#' @return index vector or match matrix
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' between(c(1,pi), 0:3,2:5+.5, T)
#' between(c(1,pi), 0:3,2:5+.5, F)
between <- function(x,
                    low,
                    high = NULL,
                    index = TRUE,
                    named = TRUE,
                    ...) {
    #browser()
    if(!is.null(dim(low)) && nrow(low == 2)){
        high <- low[2,]
        low <- low[1,]
    }
    # greater or equal to `low`
    ge.low <- outer(x, low, `>=`)
    # less or equal to `high`
    le.high <- outer(x, high, `<=`)
    inbetween <- ge.low & le.high
    if(index){
        idx.mat <- inbetween %*% 1:length(low)
        if(named){
            dimnames(idx.mat) <- list(x, 'idx')
        }
        if(any(apply(idx.mat,1, function(x)sum(x>0)) > 1)){
            # at least one x in more than one range
            ret.val <- apply(idx.mat, 1, function(v)list(v[v>0]))
        } else {
            # only unique (or no) matches
            ret.val <- apply(idx.mat, 1, sum)
        }
    } else {
        ret.val <- inbetween
        if(named){
            rownames(ret.val) <- x
            colnames(ret.val) <- paste(low, high, sep = '-')
        }
    }
    return(ret.val)
}

#' Concatenate a list of tables
#'
#' Tables are
#' The individual tables---assumed to be output from \code{table}---in
#' \code{tbl.list} are split into their name and value part which are combined
#' to 2-row matrices via \code{rbind}. Those matrices are in turn concatenated
#' via \code{cbind}.
#'
#' Names of the sub-tables are assumed to be numeric as is the case for the table
#' lists
#' @param tbl.lst list of tables
#' @param ... ignored
#'
#' @return matrix corresponding to the \code{cbind}-concatenation of \code{tbl.lst}
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' concat.tbl.list(list(table(sample(1:3, 10, T)), table(sample(7:9, 10, T))))
concat.tbl.list <- function(tbl.lst,
                            ...){
    tbl.mat <- do.call(cbind,
                       lapply(tbl.lst,
                              function(tbl){
                                  rbind(as.numeric(names(tbl)),
                                        as.vector(tbl))
                                  }
                              )
                       )
    return(tbl.mat)
}

#' vector levels
#'
#' count the number of unique entries in vector (usually an index vector with
#' integer elements)
#' @param v vector
#'
#' @return number of unique lements in \code{v}
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' v <- c(1:4, 1:3, 1:2, 1)
#' vlevels(v)
#' table(v)
vlevels <- function(v){
    return(length(unique(v)))
}

#' matrix to list
#'
#' Convert matrix \code{mat} to a list of either column vectors (for \code{which == 2})
#' or row vectors  (\code{which == 1}).
#' @param mat input matrix
#' @param which 1 for row, 2 for column vectors
#' @param ... ignored
#'
#' @return list with matrix split in vectors
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' m <- matrix(1:6, ncol = 2)
#' mat2list(m)
#' mat2list(m, 1)
mat2list <- function(mat, which = 2, ...){
    if(which == 2){
        lapply(1:ncol(mat), function(i)mat[,i])
    } else {
        lapply(1:nrow(mat), function(i)mat[i,])
    }
}

#' convert 3-element RGB-vector to color
#'
#' Wrapper for \code{\link{rgb}} that allows using it with a 3-element vector
#' \code{[r, g, b]}.
#'
#' If \code{maxColorValue} is not provided the appropriate value is guessed from
#' the values of \code{v}. If any is larger than 1, 255 is assumed. In most cases
#' that should be fine.
#'
#' @param v 3-element RGB-vector
#' @param maxColorValue passed on to \code{\link{rgb}}
#' @param ... ignored
#'
#' @return corresponding color
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' vec2rgb(1:3*80)
vec2rgb <- function(v,
                    maxColorValue = ifelse(any(v>1),255,1),
                    ...){
    if (any(v>1)) v <- v/255
    return(rgb(v[1], v[2], v[3],
               maxColorValue = maxColorValue))
}

#' Vector to HSV conversion
#'
#' Wrapper for \code{\link[grDevices]{hsv}} that allows using it with a 3-element vector
#' \code{[h, s, v]}.
#'
#' Rather than returning a hex representaion as done in \code{\link[grDevices]{hsv}}.
#' @param v a 3-element vector with the values for \code{h}, \code{s}, and \code{v}
#' @param ... ignored
#'
#' @return an RGB vector representing the color
#' @export
#'
#' @examples
#' # convert RGB colors to HSV and back
#' show.lut(apply(apply(col2rgb(rainbow(6, e=5/6)), 2, rgb2hsv), 2, vec2hsv))
vec2hsv<- function(v,
                   ...){
    if (any(v>1)) stop("values > 1 found")
    return(col2rgb(hsv(v[1], v[2], v[3])))
}

#' LAB-based Rainbow colorramp function
#'
#' Calculating the rainbow colors in LAB colorspace gives slightly better distinction
#' between neighboring colors.
#'
#' For normal usage \code{full=FALSE} should yield good results. If, however, a
#' correspondence to hue is desired,  \code{full=TRUE} should be used.
#' @param full whether to go all around to red or only to purple
#'
#' @return color ramp function
#' @importFrom grDevices colorRamp
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' rainbow_lab_ramp()(0)       # red
#' rainbow_lab_ramp()(1)       # purple
#' rainbow_lab_ramp(TRUE)(1)   # red
#' show.lut(t(rainbow_lab_ramp()(0:11/11)))
#' show.lut(t(rainbow_lab_ramp(TRUE)(0:11/11)))
rainbow_lab_ramp <- function(full = FALSE){
    base.colors <- c('red', 'yellow', 'green', 'cyan', 'blue', "#ff00ff")
    if(full) base.colors <- c(base.colors, 'red')
    return(colorRamp(base.colors,
                     space = 'Lab'))
}

#' rainbow via Lab-based colorRamp
#'
#' @param n number of colors to return
#' @param full if set, go from red to red, otherwise red to purple
#' @param ... ignored
#'
#' @return
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' rainbow_lab(12)
rainbow_lab <- function(n,
                        full = FALSE,
                        rlf  = rainbow_lab_ramp,
                        ...){
    # rainbow lab function
    return(t(rlf(seq(0,1, length.out = n))))
}

#' Split a hue range according to \code{split.tbl}
#'
#' A hue range is defined by a lower and an upper hue bound, given either as
#' 2-element vector or a \eqn{2xN} matrix where each column defines such a range.
#'
#' The splits are defined by \code{split.tbl}, a vector giving the the relative
#' weights (corresponding to "width" in colorspace) to be given to the individual
#' splits.
#'
#' The spacing between split ranges is defined by \code{blank} where a fraction
#' of \code{blank} of the original range is divided among the splits. Thus, for
#' an \eqn{n}-element vector \code{split.tbl}, each split covers \eqn{blank/(n-1)}
#' of the range and the remaining "area" is dived among the \eqn{n} splits
#' according to their weights speecified in \code{split.tbl}. See the examples.
#'
#' In either case, \code{blank} is assumed to be a "global" parameter across all
#' hue ranges.
#' @param hue.range hue range vector or matrix
#' @param split.tbl vector giving the weights of the splits (or list of such
#'                  vectors applied to corresponding matrix colunm)
#' @param blank which fraction of the range should be left blank
#'
#' @return a hue range vector (only one split) or matrix
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' hm2 <- split.hue.range(c(0,5/6), 4:2) # matrix(c(0, 4,5,11,12,20)/24, nr = 2)
#' ##           [,1]      [,2]      [,3]
#' ## [1,] 0.0000000 0.2083333 0.5000000
#' ## [2,] 0.1666667 0.4583333 0.8333333
#' hm3 <- split.hue.range(hm2, list(c(2,1,1), 3, c(1,1)))
#' ##       [,1]       [,2]      [,3]      [,4] [,5]      [,6]
#' ## [1,] 0.000 0.08333333 0.1291667 0.2083333 0.50 0.6833333
#' ## [2,] 0.075 0.12083333 0.1666667 0.4333333 0.65 0.8333333
split.hue.range <- function(hue.range,
                            split.tbl,
                            blank = 0.1) {

    if(is.null(dim(hue.range))) hue.range <- matrix(hue.range, nrow = 2)
    if(!is.list(split.tbl)) split.tbl <- list(split.tbl)
    split.tbl <- lapply(split.tbl, as.vector)    # take care of, e.g., table
    lens <- lapply(split.tbl, length)

    relative.widths <- lapply(split.tbl,                    # list
                              function(n) n/sum(n))
    diffs <- apply(hue.range, 2, diff)                      # original width
                                                            # scalar or vector
    new.ranges <- diffs * ifelse(lens==1,1,(1-blank))       # new width (- blank)
                                                            # scalar or vector
    widths <- mapply(`*`,
                     as.list(new.ranges),
                     relative.widths,
                     SIMPLIFY = FALSE)
    blanks <- mapply(function(d,n) d*blank/(length(n)-1),
                     diffs,
                     split.tbl)
    #browser()

    offsets <- mapply(function(w,b,l) c(0, cumsum(w+b))[1:l],
                      as.list(widths),
                      as.list(blanks),
                      as.list(lens),
                      SIMPLIFY = FALSE)

    ends <- mapply(`+`,
                   offsets,
                   widths,
                   SIMPLIFY = FALSE)
    rel.ranges <- mapply(rbind,
                         offsets,
                         ends,
                         SIMPLIFY = FALSE)
    real.ranges <- mapply(`+`,
                          rel.ranges,
                          as.list(hue.range[1,]),
                          SIMPLIFY = FALSE)
    if(length(real.ranges)>1){
        real.ranges <- do.call(cbind, real.ranges)
    } else {
        real.ranges <- real.ranges[[1]]
    }
    return(real.ranges)
}


#' full hue range set for cutree
#'
#' Create a set of hue ranges related through the levels of \code{cutree}.
#'
#' @param cutree result of \code{\link{cutree}}
#' @param ... passed to \code{\link{split.hue.range}}
#'
#' @return list of related hue ranges corresponding to the cut levels
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' tree.ranges(dummy.tree())
tree.ranges <- function(cutree,
                        blank = 0.1,
                        ...){
    ct.levels <- apply(cutree, 2, vlevels)
    #browser()
    ocl <- order(ct.levels)         # to order cutree by levels (if needed)
    ordered <- TRUE
    if (!identical(ocl, 1:ncol(cutree))){
        ordered <- FALSE
        rcl <- rank(ct.levels)      # to return to original order afterwards
        cutree <- cutree[, ocl]     # order columns
        ct.levels <- ct.levels[ocl]
    }
    add.top <- FALSE
    if(!any(ct.levels==1)){
        add.top <- TRUE
        cutree <- cbind(rep(1, nrow(cutree)), cutree)
    }
    # initialize with full range
    hue.splits <- list(`1`=  matrix(c( 0, 5/6),
                                    nrow = 2))
    sub.tables <- list(`1` = table(1))
    # iteratively finer
    for (i in 2:ncol(cutree)){
        # browser(text = 'level 15',
        #         expr = i == 16)
        sub.tables[[i]] <- lapply(1:vlevels(cutree[,i-1]),
                                  function(l)table(cutree[ cutree[,i-1] == l, i]))

        ## CAVEAT: the order of the sub-tables is not related to the order in
        ##         the dendrogram (arbitrary anyway) but is internally determined
        ##         by the clustering procedure.
        ##         - The index (position) corresponds to the position of the
        ##           parent cluster on the next higher level
        ##         - the level names for each sub-table refer to the corresponding
        ##           cluster indices on the current level

        # print(sub.tables[[i]])
        # print(hue.splits[[i-1]])
        #browser()
        hue.splits[[i]] <- split.hue.range(hue.splits[[i-1]],
                                           sub.tables[[i]],
                                           blank = blank^(1+.1*(i-1)),
                                           ...)[,order(concat.tbl.list(sub.tables[[i]])[1,])]
    }
    if(add.top){
        hue.splits <- hue.splits[-1]
        sub.tables <- sub.tables[-1]
    }
    names(hue.splits) <- ct.levels
    attr(hue.splits, 'sub.tables') <- sub.tables   # may need to be reordered
    class(hue.splits) <- 'tree.ranges'
    return(invisible(if(ordered) hue.splits else hue.splits[rcl]))
}

##----  Display functions ----

#' show color matrix
#'
#' show colors represented in color vector \code{cv} as matrix of colored boxes.
#'
#' @param cv color vector (actually 3xn matrix with one column per color (RGB))
#'
#' @return none
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' show.colmat(rainbow(12, end = 0.8))
show.colmat <- function(cv){
    if(is.null(dim(cv))){
        cv <- col2rgb(cv)
    }
    lcv <- ncol(cv)
    n <- ceiling(sqrt(lcv))
    m <- ifelse(lcv <= n*(n-1), n-1, n)
    plot(c(0,n)+.5, c(0, m)+.5,
         type = 'n',
         axes = FALSE,
         xlab = '',
         ylab = '')
    for (i in 1:m){             # rows
        for (j in 1:n){         # columns
            idx <- (i-1)*n + j
            cat(idx,'\n')
            if(idx <= lcv){
                print(cv[,idx])
                rect(j-.5, i-.5, j+.5, i+.5, col = vec2rgb(cv[,idx]))
            }
        }
    }
}

#' initialize hue range plot
#'
#' empty "frame" (just x-axis) with x- and y-range 0:1
#'
#' @param ... ignored
#'
#' @return none
#' @importFrom graphics axis
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' init.huerange.plot()
init.huerange.plot <- function(...){
    plot(0:1, 0:1,
         type ='n',
         axes = FALSE,
         xlab = 'hue',
         ylab = '',
         ...)
    axis(1)
}

#' convert hue range matrix to corresponding colors
#'
#' A hue range matrix is a \eqn{2\times N} matrix where each column defines a range in
#' @param hr hue range matrix
#' @param extractfun function to calculate efective hue from hue range
#' @param ... passed to \code{\link{hsv}}
#'
#' @return color (vector)
#' @importFrom grDevices hsv
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' hr <- matrix(c(0.0, 0.6,
#'                0.2, 0.7), nrow = 2, byrow = TRUE)
#' hue.range.colors(hr)
hue.range.colors <- function(hr,
                             extractfun = mean,
                             only.hues = FALSE,
                             ...){
    hues <- apply(hr, 2, extractfun)
    return(if(only.hues) hues else hsv(hues))
}


#' plot hue range
#'
#'
#' for now only horizontal
#' @param hue.range hue range or list thereof
#' @param y y-level at which to plot (scalar or vector of length of \code{hue.range})
#' @param ... passed on to \code{\link{segments}}
#' @param col consume any color that may be provided
#'
#' @return none
#' @importFrom graphics segments
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' hue.range.lines(matrix(0:5/6, nr = 2))
hue.range.lines <- function(hue.range,
                            y   = 1,
                            add = FALSE,
                            ...,
                            col = 'red'){
    ctr <- apply(hue.range, 2, mean)
    if(!add){
        init.huerange.plot(...)
    }
    hues <- apply(hue.range, 2, mean)
    # print(hues)
    segments(x0  = hue.range[1,],
             y0  = y,
             x1  = hue.range[2,],
             col = hue.range.colors(hue.range, ...),
             ...)
}

#' plot hues for a \code{\link{tree.ranges}} result
#'
#' @param tr output of \code{\link{tree.ranges}}
#' @param ... passed to \code{\link{hue.range.lines}}
#'
#' @return none
#' @importFrom graphics segments
#' @export
#'
#' @examples
#' tr <- tree.ranges(dummy.tree())
#' plot.tree.ranges(tree.ranges(dummy.tree()))
plot.tree.ranges <- function(tr,
                             add = FALSE,
                             show.tree = FALSE,
                             ...){
    n <- length(tr)
    if (! add) init.huerange.plot(...)

    line.levels <- 1 - (1:n-1)/n
    #browser()
    if(show.tree){
        for (i in 2:n){
            up.hues <- hue.range.colors(tr[[i-1]], only.hues = TRUE)
            hues <-  hue.range.colors(tr[[i]], only.hues = TRUE)
            corresp <- between(hues, tr[[i-1]])
            x.up <- up.hues[corresp]
            segments(hues, line.levels[i], x.up, line.levels[i-1],
                     col = 'lightgrey')
        }
    }
    for(i in 1:n){
        hue.range.lines(tr[[i]],
                        y   = 1 - (i-1)/n,
                        add = TRUE,
                        ...)

    }
}

#' show LUT ([color] lookup table)
#'
#' Give a visual representation of the colors represented in a LUT.
#'
#' First a color band showing the colors is presented, below the RGB components are plotted
#'
#' @param lut lookup table
#' @param colorspace
#' @param ...
#'
#' @return none
#' @export
#'
#' @examples
#' show.lut(rainbow(12, end = 11/12, v  = 0.9))
show.lut <- function(lut,
                     colorspace = c('RGB', 'HSV'),
                     ...){
    colorspace <- toupper(colorspace)
    colorspace <- match.arg(colorspace)
    if(is.null(dim(lut))) lut <- switch (colorspace,
                                         RGB = col2rgb(lut),
                                         HSV = hsv()
    )
    #print(lut)
    x <- 1:ncol(lut)
    scale <- ifelse(any(lut>1), 255,1)
    matplot(t(lut)/scale,
            type = 'l',
            col  = c('red', 'green', 'blue'),
            xlim = c(0, ncol(lut)+1),
            ylim = c(0, 1.2),
            las = 1,
            xlab = 'color index',
            ylab = colorspace,
            axes = FALSE)
    axis(1, at = x, labels = x, col = NA, col.ticks = grey(.8))
    crange <- 0:5/5
    axis(2, at = crange, labels = crange, las = 1, col = NA, col.ticks = grey(.5))
    abline(h   = crange,
           lty = 3,
           col = grey(c(.5, rep(.8, 4), .5), 0.8))
    rect(x - .5, 1.05, x + .5, 1.2,
         col = rgb(lut[1,], lut[2,], lut[3,],
                   maxColorValue = scale),
         border = NA)
}

#' Show shading in histogram-style plot
#'
#' The color matrix \code{shades} is expected to have an attribute 'reps' indicating
#' the number of shades for each main color (and, thus, the number of main colors.)
#'
#' The "base" colors are shown along the axis indicated by \code{orient}, the shades
#' are stacked perpendicular to that. If the shades were created with \code{dir = 'all'},
#' the shade stacks are aligned on the base colors
#'
#' @param shades shades as produced by \code{\link{color.shades}}
#' @param orient orientation for main colors
#' @param r "radius" for color boxes, should be less than 0.5
#' @param ... ignored
#'
#' @return none
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' rainbow.shades <- color.shades(rainbow(20, end = 0.8), sample(1:10, 20, T))
#' show.shades(rainbow.shades)
show.shades <- function(shades,
                        orient = c('horizontal', 'vertical'),
                        r      = 0.45,
                        ...){
    orient <- match.arg(orient)
    hor <- orient == 'horizontal'
    reps <- attr(shades, 'reps')
    dir <- attr(shades, 'dir')
    center <- (dir == 'all')
    if(is.null(reps)) stop('no attribute "reps"')
    s.offset <- function(n, center){
        ifelse(rep(center, length(n)),
               1 + n %/% (-2),
               0)
    }
    if(center){
        max.rep <- max(reps)
    }
    offset <- s.offset(max(reps), center)
    xmn <- ifelse(hor, 0, offset)
    ymn <- ifelse(hor, offset, 0)
    xmx <- ifelse(hor, length(reps), max(reps)+offset)
    ymx <- ifelse(hor, max(reps)+offset, length(reps))
    plot(c(xmn, xmx)+.5, c(ymn, ymx)+.5,
         type = 'n',
         axes = FALSE,
         xlab = '',
         ylab = '')
    cidx <-  1
    if(center){
        if(hor) abline(h=1, lty = 3, col =  grey(.5, 0.8)) else  abline(v=1, lty = 3, col =  grey(.5, 0.8))
    }
    for(main in 1:length(reps)){
        shade.offset <- s.offset(reps[main], center)
        for (shade in 1:reps[main]){
            xb <- ifelse(hor, main, shade+shade.offset)
            yb <- ifelse(hor, shade+shade.offset, main)
            rect(xb - r, yb - r, xb + r, yb + r,
                 col = vec2rgb(shades[, cidx]),
                 border = NA)
            cidx <- cidx + 1
        }
    }
}



#' show brain lookuptable
#'
#' show one of the lookup tables created through \code{\link{treeluts}}
#'
#' @param n number of cuts the table should have
#' @param clusters
#'
#' @return
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' \dontrun{
#'   show.brain.lut(10)
#' }
show.brain.lut <- function(n, 
                           lut = NULL,
                           clusters = read.tree()){
    if(is.null(lut)){
        # try default: file in LUTdir()
        lut.file <- sprintf(file.path(LUTdir(),
                                      'lut%03d.lut'),
                            n)
        if(file.exists(lut.file)){
            l <- readlut(lut.file, n)
        } else {
            stop("no matching file")
        }
    } else {
        if(is.matrix(lut)){
            if(nrow(lut == 3)){
                # assume directly provided LUT
                l <- lut
            } else { 
                if (ncol(lut == 3)){
                    # assume transposed LUT
                    # l <- t(lut)
                } else {
                    stop("wrong LUT format")
                }
            }
        } else {
            if (is.character(lut)){
                if(file.exists(lut)){
                    l <- readlut(lut)
                } else {
                    stop("cannot interpret '", lut, "'")
                }
            } else {
                stop("unknown type of lut")
            }
        }
    }
    # rt <- rank(clusters[,(n/5)-1])
    # rt[duplicated(rt)]<-NA
    # orrt <- order(rank(rt,
    #                    na.last = 'keep'))
    # recovered.lut <- l[orrt[1:n],]
    recovered.lut <- l[, order(reidx.cut(clusters[, (n/5)-1]))]

    lsshow <- cbind(l[,1:150], recovered.lut)
    show.lut(lsshow)
    abline(v=c(146,150))
    return(invisible(l))
}


#' show cuts
#'
#' @param cut vector with cluster assignments
#'
#' @return none
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
show.cut <- function(cut){
    n <- length(cut)
    raw <- seq_along(cut)
    n.cut <- length(unique(cut))
    plot(c(0, n+1),
         c(0, 1),
         type = 'n')
    scaled.cut <- (1:n.cut - 0.5) * n/n.cut
    points(raw, rep(0, n))
    points(scaled.cut, rep(1, n.cut))
    segments(raw, 0, scaled.cut[cut], 1)
}

## special helpers ----
#' randomize cluster assignments
#'
#' In R's implementation of \code{\link[stats]{cutree}} the first element is
#' always assigned to cluster 1, unlike what Matlab produces. This function
#' attemps to make a \code{\link[stats]{cutree}} result more similar to it's
#' Matlab counterpart (or what \code{\link{read.tree}} produces of it).
#'
#' All the columns where the number of levels does not equal the number of rows
#' (elements) is reordered (and it is considered a problem if no full column is
#' found).
#' @param cutree result of  \code{\link[stats]{cutree}}
#' @param ... ignored
#'
#' @return reordered matrix
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' ct <- dummy.tree()
#' cbind(ct,randomize.cutree(ct))
randomize.cutree <- function(cutree, ...){
    levels <- apply(cutree, 2, function(v)length(unique(v)))
    full <- levels == nrow(cutree)
    if (!any(full)){
        stop("no full cut found")
    }
    for (i in 1:length(full)){
        if (!full[i]){
            rmp <- 1:levels[i]
            # use prob to make it less likely that 1 is first
            cutree[,i] <- sample(rmp, prob = rmp)[cutree[,i]]
        }
    }
    return(cutree)
}

#' dummy tree data for examples
#'
#' The number shoud be reasonably small to keep the example results clear.
#'
#' The \code{cut} vector specifies the number of branches in the dendrogram at
#' which the cuts should be made.
#'
#' @param n number of classes [9]
#' @param cuts vector of cut numbers
#' @param ... ignored
#'
#' @return cutree result
#' @importFrom stats rnorm hclust dist cutree
#' @export
#'
#' @examples
dummy.tree <- function(n = 9,
                       cuts = round(seq(0, n,
                                        length.out = ceiling(sqrt(n+1)))[-1]),
                       ...){
    set.seed(1234)                                        # make reproducible
    dummy.data <- matrix(rnorm(n*1000),
                         nrow = n)
    hc <- hclust(dist(dummy.data))
    dt <- cutree(hc,
                 k = cuts)
    attr(dt, 'hc') <- hc
    return(dt)
}

#' show dendrogram with cut levels overlaid
#'
#' Mainly for debugging
#'
#' This needs a "real" hclust object, the output of \code{\link{cutree}} will,
#' unfortunately,  not work.
#'
#' The cut levels that are shown are set halfway between the heights of the
#'  merging levels.
#' @param hc output of \code{\link{hclust}}
#' @param cuts cut levels to show (cluster numbers)
#' @param ... passed to \code{\link{abline}} (\code{col}, \code{lty}, ...)
#'
#' @return invisible the cut levels corresponding to \code{cuts}
#' @importFrom graphics abline
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' dt <- dummy.tree()
#' dend.with.cuts(attr(dt, 'hc'), as.numeric(colnames(dt)), col = "#ffa050")
dend.with.cuts <- function(hc, cuts, ...){
    # simple dendrogram plot
    plot(hc)
    # determine appropriate cur levels
    hr <- range(hc$height)
    n <- length(hc$height)
    hrx <- diff(hr)/(2*n)
    hgt <- rev(c(min(hc$height)-hrx,    # add pseudo levels at top ...
                 hc$height,
                 max(hc$height)+hrx))   # ... and bottom
    cuts <- cuts[cuts>0 & cuts <= n+1]
    cutlvl <- (hgt[cuts+1]+hgt[cuts])/2
    # show them
    abline(h=cutlvl, ...)
    return(invisible(cutlvl))
}
