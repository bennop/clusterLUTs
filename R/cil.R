#' Colors in lookup table
#'
#' Extract the unique colors in a lookup table or palette
#' 
#' 
#' @param lut lookup table (3-row RGB color matrix [0-255])
#' @param sort 
#' @param ... passed to \code{\link{vec2rgb}} and \code{\link{rgb2hsv}} 
#'
#' @return List wih two elements
#'   - colors: unique color in LUT (3xn matrix)
#'   - cindex: index vector (length of lut) indicating for each lut entry which
#'             position in \code{colors} it corresponds to
#' @export
#'
#' @examples
colorsinlut <- function(lut,
                        sort = FALSE,
                        ...){
    if(is.matrix(lut)){
        dims <- dim(lut)
        if(!any(dims==3)){
            stop("incompatible color matrix")
        } else {
            if(dims[1] != 3){
                lut <- t(lut)
            }
        }
    }
    if(is.character(lut)){
        lut <- col2rgb(lut)
    }
    rgbs <- apply(lut, 2, vec2rgb, ...)
    n <- length(rgbs)                               # number of LUT entries
    hsvs <- rgb2hsv(lut, ...)
    uri <- !duplicated(rgbs)                        # unique RGB indices
    uhi <- !duplicated(hsvs[1,])                    # unique hue indices
    if(sum(lut[,1])==0){
        uri[1] <- uhi[1] <- FALSE                   # ignore initial BLACK entry
        if(length(rgbs)==256 && rgbs[256]=='white'){
            uhi[256] <- uri[256] <- FALSE           # ignore initial WHITE entry
        }
    }
    cat("unique hues: ", sum(uhi), '\n')
    wuhi <- which(uhi)
    unique.colors <- rgbs[wuhi]
    ordered.cols <- rgbs[wuhi[order(hsvs[1, uhi])]]
    colors <- if(sort)  ordered.cols  else  unique.colors
    color.index <- rep(0, n)
    browser()
    # show.colmat(c(rgbs, unique.colors, ordered.cols), w=n)
    
    for(i in 1:length(unique.colors)){
        color.index[which(rgbs==colors[i])] <- i
    }
    return(list(colors = col2rgb(colors),
                cindex = color.index))
}

#' Simple wrapper for \code{\link{color.shades}}
#'
#' Order of first two parameters swapped, otherwise identical to \code{\link{color.shades}}.
#' 
#' The shades are ordered dark to bright by base color in \code{inpal}
#' @param inpal color matrix or vector suitable as input to \code{\link{color.shades}}
#' @param reps vector of repetitions
#' @param ... 
#'
#' @return palette with shades added
#' @export
#'
#' @examples
#' colors <- c('red', 'green', 'blue')
#' rgbm <- col2rgb(colors)
#' repl <- c(2, 1, 1, 2, 1, 1, 3, 1, 1, 2)         # replication indices
#' rgbsm <- shade.palette(rgbm, table(repl), gscale = FALSE)
#' show.colmat(cbind(rgbm[, repl],                 # replicate orig colors
#'                  rgbsm[, rank(repl, t='r')]),   # rearrange shades
#'             width = 10)
shade.palette <- function(inpal, 
                          reps, 
                          ...){
    shade.pal <- color.shades(reps = reps,
                              col = inpal,
                              ...)
    return(shade.pal)
}

#' Reverse rank
#'
#' Similar to \code{\link[base]{rank}} but gives rank 1 to biggest entry.
#' 
#' Calculate \code{\link[base]{rank}} and subtract from \eqn{N+1} where \eqn{N} 
#' is the number of elements in \code{x}.
#' 
#' The default \code{ties.method} is required, thus, this parameter help preventing 
#' a user from using anything else. To override this \code{force} needs to be set
#' to \code{TRUE}.
#' 
#' With \code{ties.method=='random'} the ranks for equal entries in \code{x} are 
#' (pseudo) random and may vary across calls (see Examples)
#' @param x vector to rank
#' @param ties.method provide default 'random' to \code{\link[base]{rank}} to 
#'                    force unique ranks
#' @param force set to \code{TRUE} to force using \code{ties.method} other than
#'              'random'
#' @param ... passed to \code{\link[base]{rank}}
#'
#' @return ranks of elements on \code{x}
#' @export
#'
#' @examples
#' set.seed(1234)                       # (for reproducibility of example only)
#' rank(c(3,2,3,1,1,2,1,2,3,1), 
#'      ties.method = 'random')         # ->     8  6  9  3  4  7  1  5 10  2
#' rev.rank(c(3,2,3,1,1,2,1,2,3,1))     # ->     1  5  2  7  8  4  9  6  3 10
#' rev.rank(c(3,2,3,1,1,2,1,2,3,1))     # ->     2  6  3 10  8  5  7  4  1  9
rev.rank <- function(x, 
                     ties.method = 'random',
                     force = FALSE,
                     ...) {
    if(ties.method != 'random' && !force){
        ties.method <-  'random'
        warning('ties.method changed back to "random"')
    }
    return(length(x) + 1 - rank(x, 
                                ties.method = ties.method,
                                ...))
}

#' Title
#'
#' @param lut 
#' @param lut.indices 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
shade.lut <- function(lut, 
                      lut.indices,
                      ...){
    reps <- table(lut.indices)
    new.indices <- rank(lut.indices, 
                        ties.method = 'random') # only option to get unique ranks
    raw.shade.lut <- shade.palette(inpal = lut, 
                                   reps = reps,
                                   ...)
    return(raw.shade.lut[, new.indices])
}



verify.shade.lut <- function(pal.or.lut, lut.indices, ...){
    spal <- shade.lut(lut = pal.or.lut,
                      lut.indices = lut.indices,
                      ...)
    mm <- cbind(pal.or.lut[,lut.indices], spal)
    show.colmat(mm, width = ncol(spal))
    
}

hue.plot <- function(colmat,
                     hue.range = c(0, 5/6),
                     min.sat = 0.1,
                     min.val = 0.1,
                     sat.range = c(1, min.sat),
                     val.range = c(min.val, 1),
                     res = 256,
                     ...){
    
    hv.mat <- matrix(0, res, res)
    hues <- seq(hue.range[1], hue.range[2], length.out = res)
    sats <- seq(sat.range[1], sat.range[2], length.out = res/2)
    vals <- seq(val.range[1], val.range[2], length.out = res/2)
    browser()
    hv <- expand.grid(h=hues, v=vals)
    hs <- expand.grid(h=hues, s=sats[-1])
    sv <- expand.grid(s=sats, v=vals)
    bg1 <- col2rgb(hsv(h=hv$h, v=hv$v, s=1, alpha = 1))
    # show.colmat(bg1, width = 256, border = NA)
    bg2 <- col2rgb(hsv(h=hs$h, s=hs$s, v=1, alpha = 1))
    sv.plane <- col2rgb(hsv(h=1, s=sv$s, v=sv$v, alpha = 1))
    #bg <- outer(hues, vals, function(h,v) col2rgb(hsv(h=h, v=v, s=1, alpha = 1)))
    bg <- cbind(bg1, bg2)
    return(invisible(bg))
    #dim(bg) <- c(2, res, res)
    # show.colmat(bg, w=256, bor=NA)
    browser()
    png('bg.png', 
        width = 516, 
        height = 256)
    image(hues, vals, bg)
}
                     
psc <- function(file = "~/Work/4philipp/BrainLUTs/colorcentroids.txt",
                ...){
    Pscols <- read.table(file, ...)
    return(rgb2hsv(t(as.matrix(Pscols))))
}

rescale2hsv <- function(colmat,
                        ...){
    if(is.character(colmat)){
        if(file.exists(colmat)){
            colmat <- psc(file = file)
        }
    }
    col.cent <- colmat
    col.cent[2:3, ] <- 1/apply(colmat[2:3,],2, max)
}