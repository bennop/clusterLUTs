require(R.matlab)

psluts <- function(outdir = '~/Work/4philipp/BrainLUTs/luts', clusters){
    lut.length <- 256     # required by MRIcron
    apply(clusters, 2, 
          function(v){
              n <- length(unique(v))
              pal <- c(rainbow(n, end = 0.65))
              fullpal <- pal[v]
              pmat <- col2rgb(c(#'#000000',      # start with black
                                fullpal,        # colors
                                rep('#000000',  # more black
                                    lut.length - nrow(clusters) - 1), 
                                '#FFFFFF')      # end with white
                              )
              imat <- as.integer(as.vector(t(pmat)))
              #imat[imat>128] <- imat[imat>128] - 256   
              file <- path.expand(file.path(outdir, sprintf("lut%03d.lut", n)))
              cat(file,"\n")
              writelut(imat,
                       file)
              return(n)
          }
    )
}
readlut <- function(file, length = 256){
    l <- readBin(file, "integer", size = 1, signed = FALSE, n = 3*length)
    return(matrix(l, ncol = 3, byrow = F))
}

read.tree <- function(file = '~/Work/4philipp/BrainLUTs/sbm_1_145_0_result_hclust_atlas.mat'){
    ci <- readMat(file)$cluster.info
    tree <- sapply(ci, function(l)l[[1]])
    
    return(tree)
}

showlut <- function(n, clusters = read.tree()){
    l <- readlut(sprintf('/Users/puetz/Work/4philipp/BrainLUTs/luts/lut%03d.lut',
                         n))
    rt <- rank(clusters[,(n/5)-1])
    rt[duplicated(rt)]<-NA
    orrt <- order(rank(rt, na.last = 'keep')) 
    recovered.lut <- l[orrt[1:n],]
   
    lsshow <- rbind(l[1:150,], recovered.lut)

    matplot(lsshow, ty='l', col=c('red', 'green', 'blue'))
    
    return(invisible(l))
}