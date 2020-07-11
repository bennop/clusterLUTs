require(R.matlab)

##-------------------------------------------------------------------------
##
##                          C A V E A T
##
## these settings define DEFAULTS that depend on my directory structure and
## need to be adjusted to work for others
##
BLdir <- '~/Work/4philipp/BrainLUTs'
LUTdir <- file.path(BLdir, 'luts')
##
##-------------------------------------------------------------------------

## source / load Rcpp helper function
## fix path problem later
# Rcpp::sourceCpp("writelut.cpp")

#' LUT generation for hierachical clustering tree
#'
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
                     outdir   = LUTdir,
                     basename = 'lut',
                     ...){
    nc <- nrow(clusters)
    lut.length <- 256     # required by MRIcron
    last.cut <- ncol(clusters)
    apply(clusters, 2,   # for each cut
          function(v){
              n <- length(unique(v))

              #pal <- c(rainbow(n, end = 0.65))
              #fullpal <- pal[v]
              fullpal <- cut.shades(v)
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
#' @param file lookup table file
#' @param length number of color entries [256]
#' @param ... ignored
#'
#' @return color matrix
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' \dontrun{
#' lut <-  readlut('test.lut')
#' }
readlut <- function(file,
                    length = 256,
                    ...){
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
#' Each of these matrices is represented as one column in the output matrix.
#'
#' @param file Matlab hclust result file (Matlab *.mat-file)
#'
#' @return a matrix with just the cluster assignments
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' \dontrun{
#'    tree <- read.tree(path_to_treefile)
#' }
read.tree <- function(file = '~/Work/4philipp/BrainLUTs/sbm_1_145_0_result_hclust_atlas.mat'){
    ci <- readMat(file)$cluster.info
    tree <- sapply(ci, function(l)l[[1]])

    return(invisible(tree))
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
vec2hsv<- function(v,
                   ...){
    if (any(v>1)) stop("values > 1 found")
    return(col2rgb(hsv(v[1], v[2], v[3])))
}


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
            print(cv[,idx])
            if(idx <= lcv)
                rect(j-.5, i-.5, j+.5, i+.5, col = vec2rgb(cv[,idx]))
        }
    }
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
#' @param col color vector; must be a valid argument to \code{\link[grDevices]{col2rgb}}.
#' @param reps vector of repetitions
#' @param scale how much of the range to white should be covered
#'
#' @return Color matrix with on entry per column
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' color.shades(c('red', 'blue'), 3)
#' color.shades(c('red', 'green', 'blue'), 2:4)
color.shades <- function(col,
                         reps,
                         scale = 0.75){
    scale <- min(0.9, max(0.1, scale))
    len.diff <- length(col) - length(reps)
    if(len.diff > 0){
        reps <- if(length(reps == 1)){
            rep(reps, length(col))
        } else {
            c(reps, rep(1, len.diff))
        }
    }
    col.reps <- lapply(seq_along(col),
                       function(cidx){
                           # the max() contruct is neded to avoid division by
                           # zero when reps[cidx]==1
                           ramp.idxs <- 0:(reps[cidx]-1)/(max(1,
                                                              reps[cidx]-1))
                           cmat <- colorRamp(c(col[cidx],
                                               'white') )(ramp.idxs * scale)
                           return(t(cmat))
                       })
    shades <- matrix(unlist(col.reps), nrow = 3)
    attr(shades, 'reps') <- reps
    return(shades)
}


#' Show shading in histogram-style plot
#'
#' The color matrix \code{shades} is expected to have an attribute 'reps' indicating
#' the number of shades for each main color (and, thus, the number of main colors.)
#'
#' The "base" colors are shown along the axis idicate dby \code{orient}, the shades
#' are stacked perpendicular to that.
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
    if(is.null(reps)) stop('no attribute "reps"')
    xmx <- ifelse(hor, length(reps), max(reps))
    ymx <- ifelse(hor, max(reps), length(reps))
    plot(c(0, xmx)+.5, c(0, ymx)+.5,
         type = 'n',
         axes = FALSE,
         xlab = '',
         ylab = '')
    cidx <-  1
    for(main in 1:length(reps)){

        for (shade in 1:reps[main]){
            xb <- ifelse(hor, main, shade)
            yb <- ifelse(hor, shade, main)
            rect(xb - r, yb - r, xb + r, yb + r,
                 col = vec2rgb(shades[, cidx]),
                 border = NA)
            cidx <- cidx + 1
        }
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
    print(lut)
    x <- 1:ncol(lut)
    scale <- ifelse(any(lut>1), 255,1)
    matplot(t(lut)/scale,
            type = 'l',
            col  = c('red', 'green', 'blue'),
            xlim = c(0, ncol(lut)+1),
            ylim = c(0, 1.2))

    rect(x - .5, 1.05, x + .5, 1.2,
         col = rgb(lut[1,], lut[2,], lut[3,], 
                   maxColorValue = scale),
         border = NA)
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
show.brain.lut <- function(n, clusters = read.tree()){
    l <- readlut(sprintf('/Users/puetz/Work/4philipp/BrainLUTs/luts/lut%03d.lut',
                         n))


    # rt <- rank(clusters[,(n/5)-1])
    # rt[duplicated(rt)]<-NA
    # orrt <- order(rank(rt,
    #                    na.last = 'keep'))
    # recovered.lut <- l[orrt[1:n],]
    recovered.lut <- l[order(reidx.cut(clusters[, (n/5)-1]))]

    lsshow <- rbind(l[1:150,], recovered.lut)
    show.lut(lsshow)
    abline(v=c(146,150))
    return(invisible(l))
}

#' re-index cut
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

#' show cuts
#'
#' @param cut
#'
#' @return
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

#' prepare shades and reorder to match cut
#'
#' For the given cut \code{cut} analyze the internal structure of each cluster
#' and use that as basis for a shading through \code{\link{color.shades}}. Then
#' reorder the shading to arrange according to cluster members.
#'
#' @param cut vector of cluster assignments
#'
#' @return color shades reordered (see details)
#' @export
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
cut.shades <- function(cut){
    tbl <- table(cut)
    cs <- color.shades(rainbow(length(tbl), end = 0.8),
                 tbl)
    cso <- cs[, order(reidx.cut(cut))]
    attr(cso, 'reps') <- attr(cs, 'reps')
    attr(cso, 'cut') <- cut
    return(cso)
}
