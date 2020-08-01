
#' local project directory
#'
#' Used as default for \code{\link{treeluts}} and \code{\link{read.tree}},
#' simply check for existence and return input (or default).
#'
#' ~~~~~~~~~~~~~~~~~~~ C A V E A T ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#'       Adjust local defaults
#'
#' these settings define DEFAULTS that depend on my directory structure and
#' need to be adjusted to work for others
#'
#' @param basedir project directory
#' @param ... ignored
#'
#' @return path to project directory
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' projdir()
projdir <- function(basedir = '~/Work/4philipp/BrainLUTs',
                    ...){
    dir.not.found <- function(dir){
        sprintf("%d dir not found - please specify explicitly", dir)
    }
    if(!dir.exists(basedir)) stop(dir.not.found('base'))

    return(basedir)

}

#' local LUT directory
#'
#' Used as default output directory for \code{\link{treeluts}}
#'
#'
#' ~~~~~~~~~~~~~~~~~~~ C A V E A T ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#'       Adjust local defaults
#'
#' these settings define DEFAULTS that depend on my directory structure and
#' need to be adjusted to work for others
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
#' \dontrun{
#'     LUTdir()
#' }
LUTdir <- function(basedir = '.',
                   subdir = 'luts',
                   ...){
    dir.not.found <- function(dir){
        sprintf("%d dir not found - please specify explicitely", dir)
    }
    if(!dir.exists(basedir)) {
        stop(dir.not.found('base'))
    }
    LUTdir <- file.path(basedir, subdir)
    if(!dir.exists(LUTdir)) {

        success <- dir.create(LUTdir)
        ifelse(success,
               warning,
               stop)("LUT subdir not found - ",
                     ifelse(success,
                            "created",
                            "attempt to create failed"))
    }
    return(LUTdir)

}


## needed packages ----
## require(R.matlab)    # for readMat()
## require(gtools)      # for even()



##---- Rcpp helper ----------------
## load Rcpp helper function
                                        #Rcpp::loadModule("clusterLUTs", TRUE)

##---  R functions -----------------

#' LUT generation for hierarchical clustering tree
#'
#' Main function of this package.
#' \code{clusters} defines the cluster assignments per cut. The number of clusters
#' is taken from the unique values per column (need not be regular).
#'
#' @details{
#' In order to show the substructure of clusters ... (\code{\link{cutshades}})
#'
#' The LUTs are named \code{sprintf("%s%03d.lut",basename,n)} where \code{n} is
#' the number of clusters in the cut
#' }
#' @param clusters matrix with cluster assignments, one column per cut
#' @param outdir where to write the LUTs to
#' @param basename prefix for LUT names ['lut']
#' @param lut.length length of the LUT (number of entries [256])
#' @param verbose set to TRUE to see the LUT filenames
#' @param ... ignored
#'
#' @return a vector with the number of color entries per LUT generated
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib clusterLUTs
#'
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' set.seed(42)
#' treeluts(cbind(sample(1:5,20,TRUE),1:20), tempdir(), "exmpl", verbose=TRUE)
treeluts <- function(clusters   = read.tree(),
                     outdir     = LUTdir(),
                     basename   = 'lut',
                     lut.length = 256,                    # as required by MRIcron
                     verbose    = getOption('verbose'),
                     ...){
    ## browser()

    n.leaf <- nrow(clusters)
    n.cuts <- ncol(clusters)

    trc <- tree.ranges(clusters)

    ## list of color matrices
    base.lut  <-  hue.range.colors( lapply(1:n.cuts,
                                           function(i) {
                                               grDevices::col2rgb(hue.range.colors(trc[[i]])[clusters[,i]])}
                                           )
                                   )

    lut.list <- apply(clusters, 2,   # for each cut
          function(v){
              n <- vlevels(v)

                                        #pal <- c(rainbow(n, end = 0.65))
                                        #fullpal <- pal[v]
              fullpal <- cutshades(v, ...)
              pmat <- cbind(grDevices::col2rgb('#000000'),                    # initial black
                            fullpal,                                          # colors
                            col.rep('#000000',                # more black
                                    lut.length - n.leaf - 2),
                            grDevices::col2rgb(ifelse(n.leaf < lut.length-1,  # extra slot?
                                                      '#FFFFFF',              # end with white
                                                      NULL)                   # don't
                            )
              )

              imat <- as.integer(as.vector(t(pmat)))
              imat[imat>128] <- imat[imat>128] - 256
              file <- path.expand(file.path(outdir,
                                            sprintf("%s%03d.lut",
                                                    basename,
                                                    n)))
              if(verbose) cat(file,"\n")
              writelut(t(imat),
                       file)
              return(imat)
          }
          )

}

#' read lookup table
#'
#' The file will be read as binary with three (unsigned) bytes per color and
#' return a color matrix with columns corresponding to the color entries (RGB).
#'
#' If \code{length} is specified, read only that many entries, otherwise all.
#' Some consistency checks on the file size are performed.
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
                  byrow = TRUE))
}

#' read cluster tree information as obtained through Matlab
#'
#' The mat-file is expected to have one component "cluster.info" which, in turn
#' contains a list of lists. Each of the latter is a \eqn{n\times 1}{nx1} matrix with cluster
#' assignments (the number of clusters at each cut is given by the number of unique
#' entries)
#'
#' Each of these matrices (vectors) is forms one column in the output matrix as
#' would be obtained by \code{\link[stats]{cutree}}. The matrix is given column
#' names indicating the number of clusters for that cut.
#'
#' @param file Matlab hclust result file (Matlab *.mat-file)
#' @param path directory where to find file
#' @param ... ignored
#'
#' @return a matrix with just the cluster assignments
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' \dontrun{
#'    tree <- read.tree(treefile, path_to_treefile)
#' }
read.tree <- function(file = system.file('extdata/sbm_1_145_0_atlas.mat',
                                         package = 'clusterLUTs'),
                      path = '.',
                      ...){
    if(!grepl("^/", file)){      # path not absolute?
        file <- file.path(path, file)
    }
    if(!file.exists(file)){
        stop(sprintf("tree file not found\n\t'%s'",
                     file))
    }
    ci <- R.matlab::readMat(file)$cluster.info
    tree <- sapply(ci, function(l)l[[1]])
    cuts <- apply(tree, 2, vlevels)
    colnames(tree) <- cuts
    return(invisible(tree))
}

#' Create a color ramp function with shades of input
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
#' csr_red <- ColorShadeRamp('red')
#' show.colmat(t(csr_red(0:6/6)), width = 7)
ColorShadeRamp <- function(col,
                           space = c('Lab', 'rgb'),
                           ...){
    space <- match.arg(space)
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
#' shades are going \code{scale} of the way to white/black.
#'
#' If \code{col} is provided and  \code{reps} is scalar it is recycled to the
#' length of \code{col}, otherwise it is filled with ones should
#' it be shorter than \code{col}.
#'
#' The shades can be constructed either towards white, black, or, symmetrically.
#' This is determined by setting \code{direction} to 'bright', 'dark', or 'all',
#' respectively. The latter is the default.
#' When \code{gscale} is TRUE, all shade steps are the same, independent
#' of the individual repetitions, otherwise each shade range covers the range
#' given by \code{scale}.
#'
#' @param reps vector of repetitions
#' @param col color vector; must be a valid argument to \code{\link[grDevices]{col2rgb}}
#'            or a \eqn{3\times N}{3xN} RGB color matrix (one column per color).
#'            see description
#' @param scale how much of the range to black/white should be covered
#' @param gscale global shade steps?
#' @param direction which way to build the shades, see Description
#' @param col.fun function to calculate default colors
#' @param ... passed to \code{\link{ColorShadeRamp}}
#'
#' @return Color matrix with on entry per column
#' @export
#'
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' color.shades(3, c('red', 'blue'))
#' show.shades(color.shades(2:4, c('red', 'green', 'blue')))
color.shades <- function(reps,
                         col       = col.fun(length(reps)),
                         scale     = ifelse(direction == 'bright', 0.7, 0.5),
                         gscale    = TRUE,
                         direction = c('all', 'bright', 'dark'),
                         col.fun   = default.hcl,
                         ...){
    ## handle parameters
    direction <- match.arg(direction)
    scale <- min(0.9, max(0.1, scale))

    ##
    if(is.null(dim(col))){
        col <- grDevices::col2rgb(col)
    }
    len.diff <- ncol(col) - length(reps)
    if(len.diff > 0){
        reps <- if(length(reps == 1)){           # scalar:
                    rep(reps, ncol(col))         #    replicate
                } else {                         # vector:
                    c(reps, rep(1, len.diff))    #    use 1 for positions not specified
                }
    }
    max.rep <- max(reps)
    smo <- 2*max.rep %/% 2 + 1          # gives smallest odd integer >= max.rep

    col.reps <- list()

    for (cidx in 1:ncol(col)){
        ## color function (black - color - white)
        csf <- ColorShadeRamp(vec2rgb(col[,cidx]), ...)

        ## the max() contruct is needed to avoid division by
        ## zero when reps[cidx]==1
        rmp.n <- ifelse(direction == 'all' && gtools::even(reps[cidx]),
                        reps[cidx]+1 ,
                        reps[cidx])
        ramp.vals <- if(gscale){
                         if(max.rep == 1){
                             0.5
                         } else {
                             switch(direction,
                                    all = {
                                        ## print(c(reps[cidx],
                                        ## (smo + 1 - reps[cidx]) %/% 2 + 1,
                                        ## (smo + 1 + reps[cidx]) %/% 2))
 #browser()
                                       v <- seq(0.25,
                                                 0.75,
                                                 length.out = smo)[((smo + 1 - reps[cidx]) %/% 2 + 1):
                                                                   ((smo + 1 + reps[cidx]) %/% 2)]

                                    },
                                    bright = seq(0.5, 1-(1-scale)*0.5,
                                                 length.out = max.rep)[1:rmp.n],
                                    dark =  seq(0.5, (1-scale)*0.5,
                                                length.out = max.rep)[1:rmp.n]
                                    )
                         }
                     } else {
                         if(rmp.n == 1){
                             0.5
                         } else {
                             switch(direction,
                                    all = {
                                        v <- seq(0.25,
                                                 0.75,
                                                 length.out = rmp.n)
                                        if(gtools::even(reps[cidx])){
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
#' Internal rank increasing within ties so that after reordering the order of ties
#' can be preserved
#'
#' the output vector is named with (1000 * original level + position) with ties for
#' that level, see example. this can handles multiplicities of up to 999 which
#' should suffice for most practical cases.
#'
#' This is useful to reorder, e.g., a color LUT from \code{\link{color.shades}}
#' according to levels.
#'
#' @param cut cut from hierarchical clustering (through \code{link[stats]{cutree}})
#'
#' @return rank vector with names to indicate level and order within ties
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' ct <- c(4, 1, 2, 3, 1, 1, 1, 1, 4, 1, 1, 3, 1, 3)
#' reidx.cut(ct)
#' ## 4001 1001 2001 3001 1002 1003 1004 1005 4002 1006 1007 3002 1008 3003
#' ##   13    1    9   10    2    3    4    5   14    6    7   11    8   12
#' ## the first occurence of level 1 is named 1001 (see details)
#' cs <- color.shades(table(ct))
#' rocs <- cs[, reidx.cut(ct)]
reidx.cut <- function(cut){
    tbl <- table(cut)
    vals <- 1:length(unique(cut))
    vals <- as.numeric(names(tbl))
    for (i in seq_along(tbl))
        cut[cut == vals[i]] <- 1000*vals[i] + 1:tbl[i]
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
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' default.rgb()
default.rgb <- function(n = 12, ...){
    grDevices::rainbow(n, end = 5/6, alpha = 1)
}

#' HCL-based rainbow
#'
#' @param n number of colors
#' @param ... passed to rainbow_hcl
#'
#' @return rainbow colors
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' default.hcl()
default.hcl <- function(n = 12, ...){
    colorspace::rainbow_hcl(n, end = 5/6 * 360)
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
    rainbow_lab(n, ...)
}

#' prepare shades and reorder to match cut
#'
#' For the given cut \code{cut} analyze the internal structure of each cluster
#' and use that as basis for a shading through \code{\link{color.shades}}. Then
#' reorder the shading to arrange according to cluster members.
#'
#' @param cut vector of cluster assignments
#' @param col.fun function mapping a number to that many colors from a palette
#' @param ... passed to \code{\link{color.shades}}
#'
#' @return color shades reordered (see details)
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' set.seed(42);
#' (cs <- cutshades(ct <-  sample(6, 15, TRUE)))
#' op <- par(mfrow=2:1)
#' show.colmat(cs, width=15)
#' show.shades(cs)
#' par(op)
cutshades <- function(cut,
                      col.fun = default.rgb,
                      ...){
    tbl <- table(cut)
    cs <- color.shades(tbl,
                       col.fun(length(tbl)),
                       ...)
    dir <- attr(cs, 'dir')
    cso <- cs[, rank(reidx.cut(cut))]
    attr(cso, 'reps') <- attr(cs, 'reps')
    attr(cso, 'cut') <- cut
    attr(cso, 'dir') <- dir
    return(cso)
}

##----  Helper functions, mostly for internal use ----

#' binary bitwise AND
#'
#' Wrapper for \code{\link[base]{bitwAnd}}, simply providing a binary interface.
#'
#' While \code{a} or \code{b} can be any type of atomic variable, results may be
#' difficult to interpret for those other than logical or integer
#'
#' @param a boolean or numeric (integer)
#' @param b boolean or numeric (integer)
#'
#' @return bitwise AND of \code{a} and \code{b}
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' 1 %&% 3
#' pi %&% 2
#' 0.5 %&% 0.75
`%&%` <- function(a,b){
    bitwAnd(a,b)
}

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
#' @param edges which of a range's edges to include in the range
#' @param ... ignored
#'
#' @return index vector or match matrix
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' between(c(1,pi), 0:3,2:5+.5, TRUE)
#' between(c(1,pi), 0:3,2:5+.5, FALSE)
between <- function(x,
                    low,
                    high = NULL,
                    index = TRUE,
                    named = TRUE,
                    edges = NULL,
                    ...) {
    if(exists("DBG"))
        browser(expr = DBG>1)   # start browser by setting `DBG=2` on command line

    if(is.logical(high) || is.character(high)){
                                        # high was not provided and some of the other parameters were given
                                        # unnamed and, thus, assigned in a wrong order
                                        #        print(c(hasArg(high), hasArg(index), hasArg(named), hasArg(edges)))
                                        #        print(c(missing(high), missing(index), missing(named), missing(edges)))
        if(is.character(high)){
            if (is.null(edges)){
                edges <- high
                high <- NULL
            }
        } else {
            stop("apparently no `high' but unnamed parameters that take it's place - please specify")
        }
    }
    edges <- match.arg(edges, c('both', 'none', 'left', 'right'))

    if(is.null(high)){
        if(!is.null(dim(low)) && nrow(low == 2)){
                                        # low is matriz with at least 2 rows: keep 1st, assign 2nd to high
            high <- low[2,]
            low <- low[1,]
        } else {
                                        # assume pointwise matching
            high <- low
            edges <- 'both'
        }
    }

    incl.r <- edges %in% c('right', 'both')
    incl.l <- edges %in% c('left', 'both')

    low.end <- outer(x, low, ifelse(incl.l, `>=`, `>`))
    high.end <- outer(x, high, ifelse(incl.r, `<=`, `<`))
    inrange <- low.end & high.end

    if(index){
        idx.mat <- inrange %*% 1:length(low)
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
        ret.val <- inrange
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
#' Parameter \code{idx} indicates whether or not to include "index" information,
#' i.e., indicate which table an entry came from (off by default)
#'
#' Names of the sub-tables are assumed to be numeric as is the case for the table
#' lists
#' @param tbl.lst list of tables
#' @param idx whether or not to include "index" information
#' @param ... ignored
#'
#' @return matrix corresponding to the \code{cbind}-concatenation of \code{tbl.lst}
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' table.list <- list(table(sample(1:3, 10, TRUE)), table(sample(7:10, 10, TRUE)))
#' concat.tbl.list(table.list)
#' concat.tbl.list(table.list, idx = TRUE)
concat.tbl.list <- function(tbl.lst,
                            idx  =  FALSE,
                            ...){
    ## list with one matrix per tbl.lst level (cut)
    el2mat <- lapply(tbl.lst,
                     function(tbl){
                         if('table'%in% class(tbl)){
                             rbind(as.numeric(names(tbl)),
                                   as.vector(tbl))
                         } else {
                             tbl
                         }
                     })
    rows.tbl <- table(sapply(el2mat, nrow))
    if(length(rows.tbl) > 1){
        stop("incompatible inputs")
    }
    tbl.mat <- do.call(cbind,
                       el2mat)
    if(is.null(rownames(tbl.mat)) &&
       nrow(tbl.mat) == 2){
        rownames(tbl.mat) <- c('downlink',  # index of table on next (i+1) level
                               'n')         # elements in cluster
    }
    if (idx){
        index <- rep(1:length(tbl.lst),
                     times = sapply(tbl.lst,
                                    function(t.o.m){  # <- table or matrix
                                        if('table'%in% class(t.o.m)){
                                            length(t.o.m)
                                        } else {
                                            ncol(t.o.m)
                                        }
                                    } )
                     )
        tbl.mat <- rbind(index = index,     # automatically assign rowname
                         tbl.mat)
    }
    return(tbl.mat)
}


#' vector levels
#'
#' count the number of unique entries in vector (usually an index vector with
#' integer elements)
#'
#' @param v vector
#'
#' @return number of unique elements in \code{v}
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

##' convert table to matrix
##'
##' Create a matrix with names and values of the input table forming a row each.
##'
##' It is assumed that \code{tbl}'s names are all numeric
##'
##' Parameter \code{uplink} allows inclusion of additional rows below the
##' defaults rows 'idx' and 'n'. If provided as matrix, \code{uplink}
##' should have rownames.
##' @title table to matrix
##' @param tbl output of \code{\link[base]{table}}
##' @param uplink extra row(s) to be included (see description)
##' @param ... ignored
##'
##' @return a 2-row matrix representing the table.
##' @export
##' @author Benno Pütz \email{puetz@@psych.mpg.de}
##'
##' @examples
##' tbl2mat(table(1:3))
tbl2mat <- function(tbl,
                    uplink = NULL,
                    ...){
    mat <- rbind(as.numeric(names(tbl)),
                 tbl)
    rownames(mat) <- c('idx', 'n')
    colnames(mat) <- mat[1, ]

    if(!is.null(uplink)) {
        n <- ncol(mat)
        dex <- dim(uplink)
        if(is.null(dex)){                             # vector?
            if(length(uplink) == n){                  # correct length?
                mat  <- rbind(mat, uplink = uplink)
            } else {
                warning("incompatible vector 'uplink' - ignored")
            }
        } else {
            if(length(dex) == 2){                     # matrix?
                if(dim(uplink)[2] == n){              # correct ncol?
                    mat  <- rbind(mat, uplink)
                } else {
                    warning("incompatible matrix 'uplink' - ignored")
                }
            } else {
                warning("incompatible array 'uplink' - ignored")
            }
        }
    }
    return(mat)
}                                       # tbl2mat

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
#' @param v 3- or 4-element RGB-vector
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
    if (is.null(dim(v))){     # vector
        if (!(length(v) %in% 3:4)){
            stop("wrong vector length")
        } else {
            return(if(length(v) == 3){
                       grDevices::rgb(v[1], v[2], v[3],
                                      maxColorValue = maxColorValue)
                   } else {
                       grDevices::rgb(v[1], v[2], v[3], v[4],
                                      maxColorValue = maxColorValue)
                   }
                   )
        }
    } else {
        if  (!(nrow(v) %in% 3:4)){
            stop("wrong matrix dimension ([3 or 4]xN)")
        } else {
            return(if(nrow(v) == 3){
                       grDevices::rgb(v[1,], v[2,], v[3,],
                                      maxColorValue = maxColorValue)
                   } else {
                       grDevices::rgb(v[1,], v[2,], v[3,], v[4,],
                                      maxColorValue = maxColorValue)
                   }
                   )
        }
    }
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
vec2hsv <- function(v,
                    ...){
    if (any(v>1)) stop("values > 1 found")
    return(grDevices::col2rgb(hsv(v[1], v[2], v[3])))
}

#' color repetition
#'
#' Repeat color \code{col} \code{n} times as columns in a \eqn{3\times n}{3xn}
#' RGB color matrix. Analog to \code{\link[base]{rep}} which is used internally.
#'
#' @param col color (valid input to \code{\link[grDevices]{col2rgb}}
#' @param n number of repetitions (should be a positive integer value)
#'
#' @return color matrix with \code{n} columns of color \code{col}
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' col.rep('red', 3)
col.rep <- function(col, n) {
    return(matrix(rep(grDevices::col2rgb(col),                # more black
                      n),
                  nrow = 3))
}

#' LAB-based Rainbow colorramp function
#'
#' Calculating the rainbow colors in LAB color space gives slightly better distinction
#' between neighboring colors.
#'
#' For normal usage \code{full=FALSE} should yield good results (hues from red
#' to purple). If, however, a correspondence to hue is desired,
#'  \code{full=TRUE} should be used to cover
#' the full hue range from [0..1].
#' @param full whether to go all around to red or only to purple
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
#' Create a Lab-based rainbow color matrix.
#' @param n number of colors to return
#' @param rlf rainbow Lab function (\code{\link{rainbow_lab_ramp}})
#' @param ... passed to \code{rlf}
#'
#' @return color matrix
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' rainbow_lab(12)
rainbow_lab <- function(n,
                        rlf  = rainbow_lab_ramp(),
                        ...){
                                        # rainbow lab function
    return(t(rlf(seq(0,1, length.out = n), ...)))
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
#' @param split.tbl vector (possibly result of \code{\link{table}}) giving the
#'                  weights of the splits (or list of such
#'                  vectors applied to corresponding matrix colunm)
#' @param blank which fraction of the range should be left blank
#'
#' @return a hue range vector (only one split) or matrix
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' hm2 <- hue.range.split(c(0,5/6), 4:2) # matrix(c(0, 4,5,11,12,20)/24, nr = 2)
#' ##           [,1]      [,2]      [,3]
#' ## [1,] 0.0000000 0.2083333 0.5000000
#' ## [2,] 0.1666667 0.4583333 0.8333333
#' hm3 <- hue.range.split(hm2, list(c(2,1,1), 3, c(1,1)))
#' ##       [,1]       [,2]      [,3]      [,4] [,5]      [,6]
#' ## [1,] 0.000 0.08333333 0.1291667 0.2083333 0.50 0.6833333
#' ## [2,] 0.075 0.12083333 0.1666667 0.4333333 0.65 0.8333333
hue.range.split <- function(hue.range,
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
#' Create a set of hue ranges related through the levels of \code{cuts}.
#'
#' When \code{blank} is non-zero, that fraction of a hue range is left unassigned
#' when splitting for the next level. This should help for better distingushable
#' colors in the resulting hue ranges.
#'
#' @param cuts cluster assignment result of \code{\link[stats]{cutree}} (or equivalent)
#' @param blank "spacing" between hues ranges when splitting (see description)
#' @param init provide initial hues - not yet implemented
#' @param ... passed to \code{\link{hue.range.split}}
#'
#' @return list of related hue ranges corresponding to the cut levels
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' tree.ranges(dummy.tree())
tree.ranges <- function(cuts,
                        blank = 0.1,
                        init  = NULL,
                        ...){

    if (!is.null(init)){
        warning("hue initialization not yet supported - value ignored\n")
    }
    ct.levels <- apply(cuts, 2, vlevels)

    ## browser()

    ocl <- order(ct.levels)         # to order cuts by levels (if needed)
    ordered <- identical(ocl, 1:ncol(cuts)) # cuts in increasing order?
    if (!ordered){
        ordered <- FALSE
        rcl <- rank(ct.levels)      # to return to original order afterwards
        cuts <- cuts[, ocl]         # order columns
        ct.levels <- ct.levels[ocl]
    }

    add.top <- !any(ct.levels==1)
    if(add.top){                    # no top-level cut -> add one
        cuts <- cbind(rep(1, nrow(cuts)), cuts)
    }

    ## initialize with full range
    hue.splits <- list(`1`=  matrix(c( 0, 5/6),
                                    nrow = 2))
    sub.tables <- list(`1` = table(1))

    ## iteratively finer
    for (i in 2:ncol(cuts)){
        ## browser(text = 'level 15',
        ##         expr = i == 16)
        sub.tables[[i]] <- lapply(1:vlevels(cuts[,i-1]),
                                  function(l)table(cuts[ cuts[,i-1] == l, i]))

        ## CAVEAT: the order of the sub-tables is not related to the order in
        ##         the dendrogram (arbitrary anyway) but is internally determined
        ##         by the clustering procedure.
        ##         - The index (position) corresponds to the position of the
        ##           parent cluster on the next higher level
        ##         - the level names for each sub-table refer to the corresponding
        ##           cluster indices on the current level

        ## print(sub.tables[[i]])
        ## print(hue.splits[[i-1]])
        ## browser()
        hue.splits[[i]] <- hue.range.split(hue.splits[[i-1]],
                                           sub.tables[[i]],
                                           blank = blank^(1+.1*(i-1)), # shrink at ervery step
                                           ...)[,order(concat.tbl.list(sub.tables[[i]])[1,])]
    }

    if(add.top){
        hue.splits <- hue.splits[-1]
        sub.tables <- sub.tables[-1]
    }
    names(hue.splits) <- ct.levels
    names(sub.tables) <- ct.levels
    attr(hue.splits, 'sub.tables') <- sub.tables
    class(hue.splits) <- c('tree.ranges', 'list')
    return(invisible(if(ordered) hue.splits else hue.splits[rcl]))
}                                       # tree.ranges


show.tree.ranges <- function(tr,
                             ...){
    concat.tbl.list(tr)
}

##----  Display functions ----

#' show color matrix
#'
#' Show colors represented in color vector \code{cv} as matrix of colored boxes.
#'
#' By setting \code{show.idx} to TRUE, the index in the color vector is shown in
#' the center of the color tiles. Assigning a vector the same length as \code{cv},
#' the elements of that vector will be shown rather than the index.
#'
#' @param cv color vector (actually \eqn{3\times n}{3xn} matrix with one column per color (RGB))
#' @param width how many color tiles across (~\eqn{\sqrt{n}}{sqrt(n)} where \eqn{n}
#'              is the length of the color vector)
#' @param show.idx add index to display
#' @param ... passed to \code{\link[graphics]{rect}}
#'
#' @return none
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' op <- par(mfrow=c(2,1))
#' show.colmat(rainbow(12, end = 0.8))
#' show.colmat(rainbow(12, end = 0.8), width = 12)
#' par(op)
show.colmat <- function(cv,
                        width    = NULL,
                        show.idx = FALSE,
                        ...){
    if(is.null(dim(cv))){
        cv <- grDevices::col2rgb(cv)
    }
    lcv <- ncol(cv)
    n <- ifelse(is.null(width),
                ceiling(sqrt(lcv)),
                width)
    m <- ceiling(lcv/n)  # ifelse(lcv <= n*(n-1), n-1, n)
    plot(c(0,n)+.5, c(0, m)+.5,
         type = 'n',
         axes = FALSE,
         xlab = '',
         ylab = '')
    xcex  <- 1
    if(!is.logical(show.idx)){
        if(length(show.idx) != lcv){
            warning("incorrect length of 'show.idx'- ignored")
            show.idx  <- FALSE
        } else {
            if (!is.character(show.idx)){
                show.idx <- as.character(show.idx)
            }
            use.text <- show.idx
            max.width <- max(graphics::strwidth(use.text))
            if(max.width>1){
                print(max.width)
                xcex <- 1/max.width
            }
            show.idx <- TRUE
        }
    } else {
        use.text <- 1:lcv
    }
    #browser()
    d  <- 0.01                  # offset for text shadow
    for (i in 1:m){             # rows
        for (j in 1:n){         # columns
            idx <- (i-1)*n + j
            if(idx <= lcv){
                graphics::rect(j-.5, i-.5, j+.5, i+.5,
                               col = vec2rgb(cv[,idx]),
                               ...)
                if (show.idx){
                    # shadow
                    graphics::text(j+d, i-d, use.text[idx],
                                   col = 'black', cex = xcex)
                    # white text
                    graphics::text(j, i, use.text[idx],
                                   col = 'white', cex = xcex)
                }
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
    suppressWarnings(
        plot(0:1, 0:1,
             type ='n',
             axes = FALSE,
             xlab = 'hue',
             ylab = '',
             ...)
    )
    axis(1)
}

#' blank plot
#'
#' just an empty plot [0,1]x[0,1] with nothing on it (no axes, no labels)
#'
#' @param ... goes to plot
#'
#' @return none
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' init.blank.plot()
init.blank.plot <- function(...){
    suppressWarnings(
        plot(0:1, 0:1,
             type ='n',
             axes = FALSE,
             xlab = '',
             ylab = '',
             ...)
    )
}

#' convert hue range matrix to corresponding colors
#'
#' A hue range matrix is a \eqn{2\times N} matrix where each column defines a
#' range in hue space [0 .. 1] - if any hue value is found to be larger than one
#' the whole matrix is divided by 255.
#'
#' @param hr hue range matrix or a list thereof
#' @param extractfun function to calculate effective hue from hue range
#' @param only.hues when set, return (scalar) hue, otherwise convert to RGB
#'                  vector with \code{link[grDevices]{hsv}}
#' @param ... passed to \code{\link[grDevices]{hsv}}
#'
#' @return color (vector) or hue(s) depending on \code{only.hues}
#' @importFrom grDevices hsv
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' hr <- matrix(c(0.0, 0.6,
#'                0.2, 0.7), nrow = 2, byrow = TRUE)
#' hue.range.colors(hr)
#' hue.range.colors(hr, only.hues = TRUE)
#' show.colmat(col2rgb(hue.range.colors(hr)))
hue.range.colors <- function(hr,
                             extractfun = mean,
                             only.hues = FALSE,
                             ...){
    if(is.list(hr)){
        return(lapply(hr, hue.range.colors,
                      extractfun,
                      only.hues))
    } else {
        if(any(hr > 1)) {
            if(any(hr > 255)){
                stop("invalid hue value")
            } else {
                hr <- hr/255
            }
        }
        hues <- apply(hr, 2, extractfun)
        return(if(only.hues) hues else hsv(hues))
    }
}

#' Show hue range colors
#'
#' Somehere between \code{\link{show.colmat}} and \code{\link{hue.range.lines}},
#' shows the unordered colors of one in a single row of rectangles.
#'
#' When given a list of hue ranges, those lines are stacked on top of one
#' another.
#'
#' @seealso  \code{\link{hue.range.lines}} or \code{\link{tree.ranges.plot}} for
#' ordered display.
#'
#' @param hr hue range or list thereof
#' @param extractfun function to assign hue to range (\code{\link[base]{mean}})
#' @param spacing relative spacing between hue lines
#' @param y position along y-axis (mainly for single hue range)
#' @param width effective width of rectangles along axes
#' @param ... ignored
#' @param wait for internal use, should not be set by caller
#' @param v.scale relative vertical scaling to adjust rectangles' height
#' @param width.x width along x-axis
#' @param width.y width along y-axis
#'
#' @return none
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' show.hue.range(matrix(0:5/6, nr = 2))
show.hue.range <- function(hr,
                           extractfun = mean,
                           spacing    = 0.1,
                           y          = NULL,
                           width      = NULL,
                           v.scale    = 1,
                           show.idx   = FALSE,
                           ...,
                           wait       = FALSE,
                           width.x    = width,
                           width.y    = width) {
#browser(expr = wait)
    if(is.list(hr)){
        init.blank.plot()

        n <- length(hr)                 # number of levels

        if(is.null(width) && (is.null(width.x) || is.null(width.y))){
            width <- rep(1/max(sapply(hr, ncol)) * (1-spacing),
                         n)
        }
        if(is.null(width.x)) {
            width.x <- rep(width, length.out = n)
        } else {
            width.x <- rep(width.x, length.out = n)
        }
        if(is.null(width.y) || width.y == 0){
            width.y <- rep(pmax(pmin(width, 1/n), 1/n),
                           length.out = n)
        } else {
            width.y <- rep(width.y,
                           length.out = n)
        }

        ys <- 1:n/n - 0.5/n             # equidistant y-levels

        ## print(rbind(width.x, width.y, ys))

        lapply(1:n,
               function(i){
                   show.hue.range(hr[[i]],
                                  extractfun = extractfun,
                                  spacing    = spacing,
                                  y          = ys[i],
                                  width.x    = width.x[i],
                                  width.y    = width.y[i],
                                  v.scale    = v.scale,
                                  wait       = TRUE,
                                  ...)
               }
               )
    } else {
        n <- ncol(hr)
        if(!wait) {
            init.blank.plot()
            if (is.null(y)) y <-  0.5
        }
        if (is.null(width.x) || width.x == 0)
            width.x <-  1/n * (1-spacing)
        if(is.null(width.y))
            width.y <- width.x
        if(exists("DBG"))
            browser(expr = DBG>0)   # start browser by setting `DBG=1` on command line

        hues <- hue.range.colors(hr,
                                 extractfun,
                                 only.hues    = TRUE)
        xs <- 1:n/n - 0.5/n
        r.x <- width.x * 0.5
        r.y <- width.y * 0.5
        graphics::rect(xs - r.x, y - r.y * v.scale, xs + r.x, y + r.y * v.scale,
                       border = NA,
                       col    = grDevices::hsv(hues, 1, 1))
    }
    ignore <- 0                         # suppress return value
}

#' plot hue range
#'
#' for now only horizontal
#'
#' @seealso \code{\link{hue.range.colors}} for an unordered display.
#'
#' @param hue.range hue range or list thereof
#' @param y y-level at which to plot (scalar or vector of length of \code{hue.range})
#' @param add when set, add to existing plot, otherwise start a new one
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
    suppressWarnings(
        segments(x0  = hue.range[1,],
                 y0  = y,
                 x1  = hue.range[2,],
                 col = hue.range.colors(hue.range, ...),
                 ...)
    )
}

#' plot hues for a \code{\link{tree.ranges}} result
#'
#' @seealso \code{\link{hue.range.colors}} for an unordered display
#'
#' @param tr output of \code{\link{tree.ranges}}
#' @param add add to existing plot or start new one
#' @param show.tree show lines indicating the tree structure
#' @param ts.col tree segment color
#' @param ts.lty tree segment line type
#' @param ts.lwd tree segment line width
#' @param ... passed to \code{\link{hue.range.lines}}
#'
#' @return none
#' @importFrom graphics segments
#' @export
#'
#' @examples
#' tr <- tree.ranges(dummy.tree())
#' tree.ranges.plot(tr, show.tree = TRUE)
tree.ranges.plot <- function(tr,
                             add       = FALSE,
                             show.tree = FALSE,
                             ts.col    = 'lightgrey',
                             ts.lty    = 1,
                             ts.lwd    = 1,
                             ...){

    if (! add) init.huerange.plot(...)

    n <- length(tr)                       # number of tree levels
    line.levels <- rev(1 - (1:n-1)/n - 0.5/n)

    #browser()
    if(show.tree){
        # trace clusters
        for (i in 2:n){
            # hue range center(s) on previous level
            up.hues <- hue.range.colors(tr[[i-1]], only.hues = TRUE)
            # hue range center(s) on current level
            hues <-  hue.range.colors(tr[[i]], only.hues = TRUE)
            # which of the preceding range(s) are curent hues in?
            corresp <- between(hues, tr[[i-1]])
            # properly align matches on preceding level to current hues
            # (not necessarily any specific order, just correspondence)
            x.up <- up.hues[corresp]
            # show lines representing those correspondences
            segments(hues, line.levels[i], x.up, line.levels[i-1],
                     lty = ts.lty,
                     lwd = ts.lwd,
                     col = ts.col)
        }
    }

    # color bars
    for(i in 1:n){
        hue.range.lines(tr[[i]],
                        y   = line.levels[i],
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
#' @param colorspace should be either 'RGB' or 'HSV'
#' @param ... ignored
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
                                         RGB = grDevices::col2rgb(lut),
                                         HSV = hsv()
    )
    #print(lut)
    x <- 1:ncol(lut)
    scale <- ifelse(any(lut>1), 255,1)
    graphics::matplot(t(lut)/scale,
                      type = 'l',
                      col  = c('red', 'green', 'blue'),
                      xlim = c(0, ncol(lut)+1),
                      ylim = c(0, 1.2),
                      las  = 1,
                      xlab = 'color index',
                      ylab = colorspace,
                      axes = FALSE)
    axis(1, at = x, labels = x, col = NA, col.ticks = grDevices::grey(.8))
    crange <- 0:5/5
    axis(2, at     = crange,
         labels    = crange,
         las       = 1,
         col       = NA,
         col.ticks = grDevices::grey(.5))
    abline(h   = crange,
           lty = 3,
           col = grDevices::grey(c(.5, rep(.8, 4), .5), 0.8))
    graphics::rect(x - .5, 1.05, x + .5, 1.2,
                   col = grDevices::rgb(lut[1,], lut[2,], lut[3,],
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
#' rainbow.shades <- color.shades(sample(1:10, 20, TRUE))
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
        if(hor) {
            abline(h = 1, lty = 3, col =  grDevices::grey(.5, 0.8))
        } else {
            abline(v = 1, lty = 3, col =  grDevices::grey(.5, 0.8))
        }
    }
    for(main in 1:length(reps)){
        shade.offset <- s.offset(reps[main], center)
        for (shade in 1:reps[main]){
            xb <- ifelse(hor, main, shade+shade.offset)
            yb <- ifelse(hor, shade+shade.offset, main)
            graphics::rect(xb - r, yb - r, xb + r, yb + r,
                           col    = vec2rgb(shades[, cidx]),
                           border = NA)
            cidx <- cidx + 1
        }
    }
}

#' show sub tables
#'
#' In the context of this package, subtable refers to either a single table or a
#' list of tables which represents chunks of a higher level table (thus subtable).
#' They are used to trace clustering trees.
#'
#' If \code{st} is a list of lists, \code{\link{show.sub.tables}} will be called
#' recursively. For convenience, it is also possible to directly provide the output
#' of \code{\link{tree.ranges}}, where the attribute \code{sub.tables} is extracted
#' and used.
#'
#' At this point the table relation across levels is not directly shown, it has
#' to be infered from the labels
#'
#' @param st subtables or list thereof (or result of  \code{\link{tree.ranges}})
#' @param which optional subset of the subtables
#' @param ... ignored
#' @param wait for internal use, should not be set by caller
#'
#' @return string representation of sub tables
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' show.sub.tables(attr(tree.ranges(dummy.tree()), "sub.tables"))
show.sub.tables <- function(st,
                            which = NULL,
                            ...,
                            wait = FALSE) {
    if('sub.tables' %in% names(attributes(st))){
        #return(show.sub.tables(attr(st, 'sub.tables')))
       st <- attr(st, 'sub.tables')
    }

    if(is.list(st[[1]])){     # nested list of subtables
        if(is.null(which))
            which <- 1:length(st)
        # iterate over levels
        full.rep <- lapply(st[which], show.sub.tables, wait = TRUE)
        ret <-  lapply(full.rep, function(sv)paste(sv, collapse = "\n"))
        lens <- (sapply(ret, nchar)-1 ) %/% 4
        addspc <- max(lens) - lens
        offset <- sapply(addspc,
                         function(n){
                             return(substr('                                                           ',
                                    1, n))
                             }
                         )
        ret <- mapply(function(o,s)sub("\n", paste0("\n", o), s),
                      offset,
                      ret)
        lapply(rev(paste0(offset, ret)),
               function(s)cat(s, '\n\n'))
        return(invisible(ret))
    } else {
        # table level
        tbl.reps <- sapply(st,
                           function(t)utils::capture.output(t)[-1])
        lvl.strings <- apply(tbl.reps,
                             1,
                             paste,
                             collapse = '-- ')
        if(!wait){
            cat(paste(lvl.strings,
                      collapse = '\n'))
        } else {
            attr(lvl.strings, "widths") <- sapply(tbl.reps[1], nchar)
        }
        return(invisible(lvl.strings))

    }
}

#' show brain lookup table
#'
#' show one of the lookup tables created through \code{\link{treeluts}}
#'
#' The LUT can, alternatively, be given through \code{lut}: either a compatible
#' matrix (3 rows or columns) or a complete file path
#'
#' @section Testing:
#' When \code{test} is set, specific strings are returned for use with
#' \code{\link[testthat]{expect_equal}}.
#'
#' @param n number of cuts the table should have
#' @param lut alternative LUT, see description
#' @param clusters matrix with cluster assignments (1 column per cut)
#' @param verbose provide feedback?
#' @param test set to give specific return strings for tests
#' @param ... ignored
#'
#' @return LUT used (invisibly)
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' \dontrun{
#'   show.brain.lut(10)
#' }
    show.brain.lut <- function(n,
                           lut = NULL,
                           clusters = read.tree(),
                           verbose = getOption('verbose'),
                           test = FALSE,
                           ...){
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
            if(nrow(lut) == 3){
                # assume directly provided LUT
                if(verbose) cat("direct lut\n")
                if(test) return("direct lut")

                l <- lut
            } else {
                if (ncol(lut) == 3){
                    # assume transposed LUT
                    if(verbose) cat("t(lut)\n")
                    if(test)return("t(lut)")
                    l <- t(lut)
                } else {
                    stop("wrong LUT format")
                }
            }
            n <- ncol(lut)      # adjust to input
        } else {
            if (is.character(lut)){
                if(verbose) cat("lut string")
                if(file.exists(lut)){
                    l <- readlut(lut)
                    if (nrow(l) != 3){
                        stop("not an acceptable LUT file")
                    } else {
                        if (test) return("LUT from specified file")
                    }
                } else {
                    stop("cannot find '", lut, "'")
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

    ### broken
    # recovered.lut <- l[, order(reidx.cut(clusters[, (n/5)-1]))]
    #
    # lsshow <- cbind(l[,1:150], recovered.lut)
    # show.lut(lsshow)
    # abline(v=c(146,150))
    #show.lut(l)
    return(invisible(l))
}



#' show cuts
#'
#' show how levels of cut relate to individual features
#'
#' @param cut vector with cluster assignments
#'
#' @return none
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' show.cut(dummy.tree()[[1]])
show.cut <- function(cut){
    n <- length(cut)
    raw <- seq_along(cut)
    n.cut <- length(unique(cut))
    plot(c(0, n+1),
         c(0, 1),
         type = 'n')
    scaled.cut <- (1:n.cut - 0.5) * n/n.cut
    graphics::points(raw, rep(0, n))
    graphics::points(scaled.cut, rep(1, n.cut))
    graphics::segments(raw, 0, scaled.cut[cut], 1)
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
#' @param seed seed for randoom numbers to make code reproducible
#'
#' @return cutree result
#' @importFrom stats rnorm hclust dist cutree
#' @export
#'
#' @examples
#' dummy.tree()
#' dummy.tree(25)
#' dend.with.cuts(dummy.tree(25))
dummy.tree <- function(n    = 9,
                       cuts = round(seq(0, n,
                                        length.out = ceiling(sqrt(n+1)))[-1]),
                       ...,
                       seed = 1234){
    set.seed(seed)                                        # make reproducible
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
#' This needs a "real" \code{\link[stats]{hclust}} object, the output of
#' \code{\link[stats]{cutree}} will, unfortunately,  not work.
#'
#' The cut levels that are shown are set halfway between the heights of the
#' merging levels.
#'
#' When \code{cuts} are not provided, an estimate close to \eqn{\sqrt{n}{sqrt(n)}
#' where \eqn{n} is the number of leaves.
#' @param dt output of \code{\link{dummy.tree}}
#' @param cuts cut levels to show (cluster numbers) (optional, see description)
#' @param cut.col color to use for cut lines ("#ffa050" which is a light orange)
#' @param ... passed to \code{\link{abline}} (\code{col}, \code{lty}, ...)
#'
#' @return invisible the cut levels corresponding to \code{cuts}
#' @importFrom graphics abline par
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' dt <- dummy.tree()
#' dend.with.cuts(dt)
dend.with.cuts <- function(dt,
                           cuts = NULL,
                           cut.col =  "#ffa050",
                           ...){
    # simple dendrogram plot
    if('hc' %in% names(attributes(dt))){
        # interpret as output of dummy.tree and extract underlying hclust object
        hc <- attr(dt, 'hc')
    } else {
        stop("not recognized as dummy.tree output")
    }

    hr <- range(hc$height)
    n <- length(hc$height)
    hrx <- diff(hr)/(2*n)
    h.low <- min(hr) - 3*hrx
    dg <- stats::as.dendrogram(hc)
    #browser()
    op <- par(mar=c(4,5,4,5)+.1)
    DEX <- FALSE
    if(!requireNamespace("dendextend")){       # standard plot
        plot(hc,
             las  = 1,
             xlab = 'cluster index',
             axes = FALSE,
             ...)
    } else {                                 # can use raise.dendrogram()
        DEX <- TRUE
        plot(dendextend::raise.dendrogram(dg, -h.low),
             axes = FALSE,
             ...)
    }
    par(op)

    if(is.null(cuts)){
        cuts <- as.numeric(colnames(dt))
    }

    # determine appropriate cur levels
    hgt <- rev(c(min(hc$height)-hrx,    # add pseudo levels at top ...
                 hc$height,
                 max(hc$height)+hrx))   # ... and bottom
    cuts <- cuts[cuts > 0 & cuts <= n+1]
    raise <-  ifelse(DEX, h.low, 0)
    cutlvl <- (hgt[cuts+1] + hgt[cuts])/2 - raise
    # show them
    abline(h = cutlvl, col = cut.col, ...)
    axis(4,
         at = cutlvl,
         labels = sub(' ', '  ',sprintf("%2d",cuts)),
         col = NA,
         lwd = 0,
         col.ticks = cut.col,
         las = 1,
         ...)
    return(invisible(cutlvl+raise))
}
