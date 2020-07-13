##----------- C A V E A T --------------------------------------------------------------
##
##      Adjust local defaults
##      
## these settings define DEFAULTS that depend on my directory structure and
## need to be adjusted to work for others
##
BLdir <- '~/Work/4philipp/BrainLUTs'
LUTdir <- file.path(BLdir, 'luts')
##
##---

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
#' @author Benno Pütz \email{puetz@@psych.mpg.de}
#'
#' @examples
#' \dontrun{
#'    tree <- read.tree(path_to_treefile)
#' }
read.tree <- function(file = '~/Work/4philipp/BrainLUTs/sbm_1_145_0_result_hclust_atlas.mat'){
    ci <- readMat(file)$cluster.info
    tree <- sapply(ci, function(l)l[[1]])
    cuts <- apply(tree, 2, function(v)length(unique(v)) )
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
#' @export
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
#' @export
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

default.rgb <- function(n, ...){
    rainbow(n, end = 5/6)
}

default.hcl <- function(n, ...){
    rainbow_hcl(n, end = 5/6)
}

default.lab <- function(n, ...){
    rainbow_lab(n, end = 5/6)
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
rainbow_lab <- function(n,
                        full = FALSE,
                        ...){
    # rainbow lab function
    return(t(rlf(seq(0,1, length.out = n))))
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
show.brain.lut <- function(n, clusters = read.tree()){
    l <- readlut(sprintf('/Users/puetz/Work/4philipp/BrainLUTs/luts/lut%03d.lut',
                         n))
    
    
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
#' m <- matrix(rnorm(15*1000), nrow = 15)
#' hc <- hclust(dist(m))
#' ct <- cutree(hc, seq(5,15, by=5))
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
    