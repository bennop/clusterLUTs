## class helpers ----
check_huerange <- function(object) {
    errors <- character()
    
    if (object@from < 0) {
        msg <- "'from' is less than 0"
        errors <- c(errors, msg)
    }
    
    if (object@from > 1) {
        msg <- "'from' is larger than 1"
        errors <- c(errors, msg)
    }
    
    if (object@to < 0) {
        msg <- "'to' is less than 0"
        errors <- c(errors, msg)
    }
    
    if (object@to > 1) {
        msg <- "'to' is larger than 1"
        errors <- c(errors, msg)
    }
    
    if (object@from > object@to) {
        msg <- paste("'from' is less than 'to' - should be ordered")
        errors <- c(errors, msg)
    }
    
    if (length(errors) == 0) TRUE else errors
}

## Class definition ----

huerange <- setClass("huerange",
                     representation(from='numeric', to='numeric'),
                     prototype(from = 0, to = 1),
                     validity = check_huerange)

## Methods ----
setMethod("print", 
          signature(x = "huerange"), 
          function(x) {
              cat(sprintf("%f - %f\n",x@from, x@to))
          })

setMethod("plot", 
          signature(x = "huerange"), 
          # simply consume any color that might be provided
          function(x, y=1, col = 'red', ...) {
              cat("hue\n")
              lines(x = c(x@from, x@to),
                    y = rep(y,2), 
                    col = hsv(center(x)),
                    ...)
          })

setGeneric("center", function(object) {
    standardGeneric("center")
})
setMethod("center", 
          signature(object = "huerange"), 
          function(object) {
              mean(c(object@from, object@to))
          })

setGeneric("split", function(x, n, blank = 0.1) {
    standardGeneric("split")
})

setMethod("split",
          signature(x = "huerange"),
          function(x, n, blank = 0.1) {
              n <- as.vector(n)    # take care of, e.g., table
              if(length(n)==1) return(x)
              relative.widths <- n/sum(n)
              new.range <- (x@to-x@from) * (1-blank)
              widths <- new.range * relative.widths
              blanks <- (x@to-x@from) * blank/(length(n)-1)
              sub.ranges <- list()
              start <- x@from
              for (i in 1:length(n)){
                  end <- start + widths[i]
                  sub.ranges <- c(sub.ranges,
                                  new('huerange',
                                      from = start,
                                      to = end))
                  start <- end + blanks
              }
              return(sub.ranges)
          })


