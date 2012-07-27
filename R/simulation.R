

.geno.equi <- function(Np, map) {
    # FIXME: what about marker bins ? one sample per bin or per marker ?
    # One per marker it is.
    ret <- sapply(1:Np, function(i) sample(c('a', 'b'), nrow(map), replace=T))
    rownames(ret) <- map$markers
    ret
}


.mosaic2.ind <- function(Ng, Lc, map, gen) {
    .debug.out("DEBUG", Ng, Lc, "\n")
    .debug.out(dim(map))
    .debug.out(dim(gen))
    Nb <- rpois(1, Ng)
    .debug.out(Nb)
    .debug.out(all(as.character(map$markers) %in% rownames(gen)))
    G <- gen[as.character(map$markers), ]
    .debug.out("Have", nrow(G), "markers.\n")
    if(Nb == 0) {
        return(G[, sample.int(ncol(G), 1)])
    }
    marker.segments <- matrix(F, ncol=nrow(map), nrow=nrow(map))
    break.pos <- sort(runif(Nb, 0, Lc))
    m <- 1
    b <- 1
    rowb <- 1
    .debug.out("Generating", Nb, "breaks.\n")
    while(b <= Nb) {
        loc.min <- map$locus[m]
        while(b < Nb && break.pos[b] < loc.min) {
            b <- b + 1
        }
        #.start <- m
        while(map$locus[m] <= break.pos[b]) {
            marker.segments[rowb, m] <- T
            m <- m + 1
        }
        #rowb <- rowb + (.start != m)
        rowb <- rowb + 1
        b <- b + 1
    }
    if(m <= nrow(map)) {
        marker.segments[rowb, m:nrow(map)] <- T
    } else {
        rowb <- rowb - 1
    }
    .debug.out("Kept", rowb, "non-empty segments.\n")
    #marker.segments <- marker.segments[1:rowb, ]
    ret <- unlist(lapply(1:rowb,
                         function(r) G[marker.segments[r,], sample.int(ncol(G), 1)]))
    .debug.out("Generated", length(ret), "markers.\n")
    ret
}


.geno.mosaic2 <- function(Np, Ng, map, gen) {
    Lc <- map[nrow(map), 'locus']
    Ng <- Lc * Ng * .01         # Lc in Morgan, not centiMorgan. 
    ret <- sapply(1:Np, function(i) .mosaic2.ind(Ng, Lc, map, gen))
    rownames(ret) <- map$markers
    ret
}


#> k <- .geno.equi(20, clu$get.haplotypes.map('chr01'))
#> h <- .geno.mosaic(50, 100, clu$get.haplotypes.map('chr01'), k)
#> j <- .geno.mosaic(50, 1000, clu$get.haplotypes.map('chr01'), k)


.mosaic.make.breaks <- function(Lc, Ng) {
    Nb <- rpois(1, Lc * Ng * 0.01)        # Lc in Morgan, not centiMorgan.
    break.pos <- sort(runif(Nb, 0, Lc))
    segment.start <- c(0, break.pos)      # break loci delimit intervals
    segment.end <- c(break.pos, Inf)      # 
    list(start=segment.start, end=segment.end)
}


.mosaic.make.segment <- function(map, gen, start.pos, end.pos, p.num) {
    markers.ok <- map[, 'locus'] >= start.pos & map[, 'locus'] < end.pos
    #print(c(start.pos, end.pos))
    #print(markers.ok)
    markers <- as.character(map$markers[markers.ok])
    #print(c(p.num, markers))
    if (length(markers) == 0) {
        return(NULL)
    }
    gen[markers, p.num]
    #rep(p.num, length(markers))
}


.mosaic.make.geno <- function(Ng, map, gen) {
    Lc <- map[nrow(map), 'locus']

    segment.list <- .mosaic.make.breaks(Lc, Ng)
    segment.list$p.num <- sample.int(ncol(gen), length(segment.list$start), replace=T)
    .debug.out(as.data.frame(segment.list))

    unlist(lapply(1:length(segment.list$p.num),
                  function(i) {
                      .mosaic.make.segment(map, gen, segment.list$start[i], segment.list$end[i], segment.list$p.num[i])
                  }))
}


.geno.mosaic <- function(Np, Ng, map, gen) {
    ret <- sapply(1:Np, function(i) .mosaic.make.geno(Ng, map, gen))
    rownames(ret) <- map$markers
    ret
}



.gen.geno.equi <- function(Np) {
    function(map, gen) .geno.equi(Np, map)
}

.gen.geno.mosaic <- function(Np, Ng) {
    function(map, gen) .geno.mosaic2(Np, Ng, map, gen)
}




.make.ancestors.genotype <- function(Na, Nm) {
    ret <- t(matrix(letters[1:Na], nrow=Na, ncol=Nm))
    rownames(ret) <- paste('M', 1:Nm, sep="")
    ret
}



.make.marker.dists <- function(Nm, len=4, shapes=c(.5), peaks=c(2), peak.weights=c(1)) {
    wmaker <- function(p, w, s) {
        f <- function(x) w*exp(-abs(x-p)^s)
        attr(f, 'p') <- p
        attr(f, 'w') <- w
        attr(f, 's') <- s
        f
    }
    if(length(peaks) == 0) {
        ret <- sort(runif(Nm, 0, len))
    } else {
        wfuns <- lapply(1:length(peaks),
                        function(i) wmaker(peaks[i], peak.weights[i], shapes[i]))
        wfun <- function(x) {
            sapply(x, function(x) sum(sapply(wfuns, function(f) f(x))))
        }
        norm.const <- integrate(wfun, 0, len)$value / Nm
        #wfun.dist <- function(x) norm.const / wfun(x)
        mrk.pos <- 0:(Nm-1)*len/Nm
        norm.const <- max(wfun(mrk.pos))
        wfun.dist <- function(x) (norm.const - wfun(x)) / norm.const

        ret <- wfun.dist(mrk.pos)
    }
    cs <- cumsum(ret)
    #cat("   len=", len, " Nm=", Nm, " p=", length(peaks), "   length=", cs[length(cs)], "\n", sep="")
    ret <- ret * 100 * len / cs[length(cs)]
    cs <- cumsum(ret)
    cat("   len=", len, " Nm=", Nm, " p=", length(peaks), "   length=", cs[length(cs)], "\n", sep="")
    ret
}


.generate.map <- function(Nc=10, Nm=rep(200, Nc)) {
    cat(Nc, "chromosomes to generate.\n")
    cat(Nm, "markers in each chromosome.\n")
    marker.names <- paste('M', 1:sum(Nm), sep='')
    marker.chrom <- unlist(lapply(1:Nc, function(i) rep(i, Nm[i])))
    chroms <- lapply(1:Nc,
                     function(i) {
                         m <- Nm[i]
                         mnames <- marker.names[marker.chrom == i]
                         cat("Generating map for chromosome", i, "with", m, "markers.\n")
                         npeaks <- rpois(1, 1)
                         len <- rnorm(1, 3, 1)
                         if(npeaks > 0) {
                             peaks <- runif(npeaks, 0, len)
                             pweights <- runif(npeaks, 0, 1)
                             pweights <- pweights / max(pweights)
                             pshapes <- 10^runif(npeaks, -1, 1)
                         } else {
                             peaks <- c()
                             pweights <- c()
                             pshapes <- c()
                         }
                         dists <- .make.marker.dists(m - 1,
                                                    len,
                                                    pshapes,
                                                    peaks,
                                                    pweights)
                         cvec <- rep('0', 1 + 2 * m)
                         cvec[1] <- paste('*chr', i, sep='')
                         cvec[2] <- m
                         cvec[1 + 1:m * 2] <- mnames
                         cvec[2 + 1:(m-1) * 2] <- dists
                         cvec
                     })
    ret <- matrix('', ncol=max(sapply(chroms, length)), nrow=length(chroms))
    for(i in 1:Nc) {
        ret[i, ] <- NA
        ret[i, 1:length(chroms[[i]])] <- as.character(chroms[[i]])
    }
    list(map=ret, marker.names=marker.names, marker.chrom=marker.chrom)
}



generate.data <- function(Nc, Nancestors, Nhaplotypes, Nm=80+sample.int(40, Nc, replace=T), Nmissing=0, mosaic.Ng=50, filename.prefix="generated") {
    map <- .generate.map(Nc, Nm)
    ancestors <- .make.ancestors.genotype(Nancestors, sum(Nm))
    chromosomes <- apply(map$map, 1, .read.chrom, 1)
    cat("Generating parental genotypes...\n")
    mosaic <- .gen.geno.mosaic(Nhaplotypes, mosaic.Ng)
    parent.genos <- lapply(1:length(chromosomes),
                           function(i) {
                               cat('... for chromosome', i, '\n')
                               mosaic(chromosomes[[i]]$map, ancestors)
                           })
    parent.genos.names <- c('Marker', paste("Parent", 1:Nhaplotypes, sep=""))
    d <- Sys.Date()
    map.filename <- paste(filename.prefix, '.', d, '.map.txt', sep='')
    geno.filename <- paste(filename.prefix, '.', d, '.geno.txt', sep='')
    log.filename <- paste(filename.prefix, '.', d, '.log', sep='')
    cat("Writing genetic maps to", map.filename, "\n")
    writeLines(apply(map$map, 1, function(r) paste(na.omit(r), collapse=" ")), map.filename, sep="\n")
    cat("Writing genotypes to", geno.filename, "\n")
    lapply(1:length(parent.genos),
           function(g) {
               m <- cbind(rownames(parent.genos[[g]]), parent.genos[[g]])
               colnames(m) <- parent.genos.names
               write.table(m,
                           geno.filename, quote=F, na='-',
                           append=g!=1, col.names=g==1, row.names=F)
           })
    parent.genos <- read.table(geno.filename, )
    rownames(parent.genos) <- parent.genos[, 1]
    parent.genos <- as.matrix(parent.genos[, -1])
    if(Nmissing > 0) {
        parent.genos[sample.int(length(parent.genos), Nmissing)] <- NA
    }
    write.table(cbind(rownames(parent.genos), parent.genos), geno.filename,
                quote=F, na="-", append=F, col.names=F, row.names=F)
    .cat <- function(...) cat(..., file=log.filename, append=T)
    .cat("Data generated on ", date(), "\n\n", sep="")
    
    .cat(Nancestors, "ancestors\n")
    .cat(Nc, "chromosomes were generated, each containing", Nm, "markers.\n")
    .cat(Nhaplotypes, " mosaic descendants with Ng=", mosaic.Ng, "\n", sep="")
    .cat(Nmissing, "missing data were randomly assigned in genotypes.\n")

    list(map=map$map, geno=parent.genos)
}
