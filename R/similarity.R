.test.loci <-function(chrom, step.size) {
    #.debug.out("dim(chrom) =", dim(chrom), "\n")
    steps <- lapply(1:nrow(chrom),
                    function(i) {
                        pos <- chrom[i, 'locus']
                        d <- chrom[i, 'dist.from']
                        if (is.infinite(d)) {
                            pos
                        } else if (d > step.size) {
                            ds <- d / step.size
                            if (ds == floor(ds)) {
                                seq(pos, pos + d - step.size, step.size)
                            } else {
                                seq(pos, pos + d, step.size)
                            }
                        } else if (d > 0) {
                            pos
                        }
                    })
    #unique(unlist(steps))
    unlist(steps)
}


.quick.find.segment <- function(seg, i, x) {
    nearest <- seg[which(x == min(x) & i == 1)]
    if (length(nearest) > 0) {
        lapply(nearest, function(s) x[seg == s & x != 0])
    }
}


.prepare.similarity.data <- function(desc.chrom, full.chrom,
                                     step.size, window.length) {
    .info.out("Computing loci for", nrow(desc.chrom), "with a step size of",
             step.size, "cM.\n")
    loci <- .test.loci(desc.chrom, step.size)
    map <- full.chrom
    double.markers.removed <- which(map[, 'dist.to'] <= MARKER.DOUBLE.EPSILON)
    if (length(double.markers.removed) > 0) {
        #map <- map[-double.markers.removed, ]
        map <- .remove.markers(map, double.markers.removed)
    }
    hws <- window.length * .5
    windows <- lapply(loci,
                      function(p) which(abs(p - map[, 'locus']) <= hws))
    x <- lapply(1:length(loci),
                function(i) abs(map[windows[[i]], 'locus'] - loci[i]))
    list(loci=loci,
         select=windows,
         x=x,
         map=map,
         density=sapply(windows, length),
         double.markers.removed=double.markers.removed)
}


.compute.similarity <- function(full.chrom, simi.data, gen, i, j, w1, w2, na.replace=F) {
    #cat("dim(gen) = [", paste(dim(gen), collapse=', '), '] i=', i, ' j=', j, '\n', sep='')
    iab <- .compute.I(gen, i, j, full.chrom, na.replace)
    #if (length(simi.data$double.markers.removed)) {
        #iab <- iab[-simi.data$double.markers.removed, ]
        #print(iab)
    #}
    M <- simi.data$map
    S <- simi.data$select
    X <- simi.data$x
    P <- simi.data$loci
    simi <- sapply(1:length(P),
                   function(n) {
                       select <- S[[n]]
                       if (length(select) == 0) { return(0) }
                       x <- X[[n]]
                       i <- iab$I[select]
                       #print(data.frame(row.names=full.chrom$markers[select], x=x, i=i))
                       i.ok <- !is.na(i)
                       if (!all(i.ok)) {
                           if (!any(i.ok)) { return(0) }
                           i <- i[i.ok]
                           select <- select[i.ok]
                           x <- x[i.ok]
                       }
                       SW1 <- sum(w1(x) * i)
                       all.w2 <- sapply(.quick.find.segment(iab$seg[select], i, x),
                                        function(xlist) sum(w2(xlist)))
                       SW2 <- ifelse(length(all.w2) > 0, max(all.w2), 0)
                       #.debug.out(P[n], "|", w1(x) * i, "|", SW2, "\n")
                       SW1 + SW2
                   })
    list(locus=simi.data$loci,
         similarity=simi,
         density=simi.data$density)
}
