
MAP.SIZE.TOLERANCE <- .1

MARKER.DOUBLE.EPSILON <- 1.e-10


.compute.dists <- function(posvec) {
    dists <- diff(posvec)
    list(to=c(Inf, dists), from=c(dists, Inf))
}


.find.doubles <- function(map) {
    starts <- which(abs(map$dist.from) <= MARKER.DOUBLE.EPSILON & abs(map$dist.to) > MARKER.DOUBLE.EPSILON)
    ends <- which(abs(map$dist.to) <= MARKER.DOUBLE.EPSILON & abs(map$dist.from) > MARKER.DOUBLE.EPSILON)
    if(length(starts) != length(ends)) {
        cat(paste("length(starts) =", length(starts), "length(ends) =", length(ends), "\n"))
        print(map)
        stop("BUG")
    }
    cbind(starts, ends)
}


.remove.markers <- function(map, markers) {
    map <- map[!as.character(map$markers) %in% markers, ]
    tf <- .compute.dists(map$locus)
    map$dist.to <- tf$to
    map$dist.from <- tf$from
    attr(map, 'doubles') <- .find.doubles(map)
    map
}


.init.map <- function(list.of.maps.as.list) {
    for (n in names(list.of.maps.as.list)) {
        map.as.list <- list.of.maps.as.list[[n]]
        dists <- .compute.dists(map.as.list$locus)
        map.as.list$dist.to <- dists$to
        map.as.list$dist.from <- dists$from
        doubles <- .find.doubles(map.as.list)
        attr(map.as.list, 'doubles') <- doubles
        attr(map.as.list, 'chrom.name') <- n
        list.of.maps.as.list[[n]] <- map.as.list
    }
    list.of.maps.as.list
}


# extracts the data for ONE chromosome from a row possibly padded with NA's.
.read.chrom <- function(chromvec, i) {
    chrname <- substring(chromvec[i], 2)   # strip leading *
    #chromvec <- chromvec[-1]   # strip chromosome name from vector
    nmark <- as.numeric(chromvec[i + 1])
    cat("reading", chrname, "/", nmark, "markers\n")
    ofs = i +  seq(1, 2 * nmark, 2)
    dists <- as.numeric(chromvec[ofs])
    mrkname <- as.vector(chromvec[ofs + 1])
    total <- dists[1]  # total number of markers
    dists <- dists[-1]  # strip the number of markers (first value), replace with 0
    dist.to <- c(Inf, dists)
    dist.from <- c(dists, Inf)
    d <- data.frame(markers=mrkname, dist.to=dist.to, dist.from=dist.from,  locus=cumsum(c(0, dists)))

    list(name=chrname, map=d, .next=i + 1 + 2 * nmark)
}


.read.map <- function(filename) {
    #dat <- as.matrix(read.table(filename, fill=NA, colClasses='character'))
    dat <- scan(filename, fill=NA, comment.char='', what='character', quiet=T)
    ret = list()
    #for (i in 1:dim(dat)[1]) {
    i <- 1
    while(i < length(dat)) {
        ch <- .read.chrom(dat, i)
        i <- ch$.next
        doubles <- .find.doubles(ch$map)
        .debug.out("Chromosome", ch$name, "has", dim(ch$map)[1], "markers\n")
        #ret[[ch$name]] <- list(name=ch$name, map=ch$map, doubles=doubles)
        map <- ch$map
        attr(map, 'chrom.name') <- ch$name
        attr(map, 'doubles') <- doubles
        ret[[ch$name]] <- map
    }
    ret
}



.read.gen <- function(filename, na.strings=c("-", "NA", ".")) {
    gen <- read.table(filename, header=T, na.strings=na.strings, colClasses='character')
    rownames(gen) <- gen[, 1]
    as.matrix(gen[, -1])
}


.compute.I <- function(gen, i, j, chrom, na.replace=F) {
    # na.replace : policy regarding NA values in genotypes :
    #       = NA : keep NA's in the table of I, and ignore NA segments when looking for nearest TRUE segment
    #       = T  : replace NA's with TRUE in the table of I
    #       = F  : replace NA's with FALSE in the table of I
    if (i==j) {
        i.table <- rep(1, nrow(chrom))
        names(i.table) <- chrom$markers
        seg <- rep(1, nrow(chrom))
        #ret <- data.frame(I=i.table, seg=seg)
        #ret <- matrix(c(i.table, seg), ncol=2)
        ret <- list(I=i.table, seg=seg)
        #attr(ret, 'chrom.name') <- chrom$name
        attr(ret, 'chrom.name') <- attr(chrom, 'chrom.name')
        attr(ret, 'i') <- i
        attr(ret, 'j') <- j
        attr(ret, 'na.replace') <- na.replace
        return(ret)
    }
    map <- chrom
    map.markers <- as.character(map$markers)
    mark <- map.markers[map.markers %in% rownames(gen)]
    #gen <- gen[as.character(map$markers), ]
    gen <- gen[mark, ]
    doubles <- attr(chrom, 'doubles')
    i.table <- gen[, i] == gen[, j]
    names(i.table) <- rownames(gen)
    # double markers :
    # inside each group of double markers, combine the values of I :
    # if all values are NA, result is NA
    # otherwise, NA's are ignored.
    # any FALSE results in all being set to FALSE (NA's included).
    # otherwise, all are set to TRUE (NA's included).
    if (!is.null(doubles) && nrow(doubles) > 0) {
        reord.doubles <- cumsum(map.markers %in% mark)
        for (d in 1:nrow(doubles)) {
            start.i <- reord.doubles[doubles[d, 1]]
            end.i <- reord.doubles[doubles[d, 2]]
            slice <- i.table[start.i:end.i]
            # if every value in slice is NA, resulting common value is NA,
            # so don't touch anything.
            if (!all(is.na(slice))) {
                i.table[start.i:end.i] <- all(slice, na.rm=T)
            }
        }
    }
    # hack NA's to perform correct segmentation
    if(!is.na(na.replace)) {
        i.table[is.na(i.table)] <- as.numeric(na.replace)
        seg = cumsum( c(0, i.table) != c(i.table, 0))[-length(i.table)-1]
    } else {
        which.na <- is.na(i.table)
        i.table[which.na] <- -Inf  # need comparable yet out-of-band value
        seg = cumsum( c(0, i.table) != c(i.table, 0))[-length(i.table)-1]
        # restore NA's
        i.table[which.na] <- NA
    }

    if (is.null(i.table)) {
        stop("BUG")
    }

    #ret <- data.frame(I=i.table, seg=seg)
    #ret <- matrix(c(i.table, seg), ncol=2)
    ret <- list(I=i.table, seg=seg)
    attr(ret, 'chrom.name') <- attr(chrom, 'chrom.name')
    attr(ret, 'i') <- i
    attr(ret, 'j') <- j
    attr(ret, 'na.replace') <- na.replace
    ret
}


#.remove.double.markers <- function(iab, map) {
#    # filter out double markers, keeping only the first in each group
#    .debug.out("dim(I) =", dim(iab), "dim(map) =", dim(map), "\n")
#
#    to.be.removed <- which(map$dist.to <= MARKER.DOUBLE.EPSILON)
#    if(length(to.be.removed) > 0) {
#       #.debug.out("Filtering out", length(to.be.removed), "double markers.\n")
#       removed <- map[to.be.removed, 'markers']
#       iab <- iab[-to.be.removed, ]
#       map <- map[-to.be.removed, ]
#       # just in case the final marker in map was removed, set final dist.from to Inf
#       map[nrow(map), 'dist.from'] <- Inf
#       #iab$seg <- iab$seg[-to.be.removed]
#    } else {
#        removed <- NULL
#    }
#    list(I=iab, map=map, removed=removed)
#}


.realign.loci <- function(ref.map, map) {
    .debug.out(head(ref.map))
    .debug.out(head(map))
    ref.markers <- ref.map$markers[ref.map$markers %in% map$markers]
    markers <- map$markers[map$markers %in% ref.map$markers]
    common.markers <- intersect(ref.markers, markers)
    .debug.out("maps have", length(common.markers), "markers in common.\n")
    submap <- map[map$markers %in% common.markers, ]
    subref.map <- ref.map[ref.map$markers %in% common.markers, ]
    if (max(abs(subref.map$locus - submap$locus)) > MAP.SIZE.TOLERANCE) {
        stop("The difference in size between the haplotypes map and MCQTL consensus genetic maps MUST NOT exceed ", MAP.SIZE.TOLERANCE, " cM.")
    }
    approx(y=subref.map$locus, x=submap$locus, xout=map$locus, rule=2)$y
}


.create.test.maps <- function(haplotypes.map, mcqtl.consensus.map) {
    haplotypes.chroms <- names(haplotypes.map)
    consensus.chroms <- names(mcqtl.consensus.map)
    chroms.to.merge <- intersect(haplotypes.chroms, consensus.chroms)
    merged.maps <- list()
    for(name in chroms.to.merge) {
        # dist.from and dist.to are meaningless after merge, so we don't keep them.
        dmap <- mcqtl.consensus.map[[name]][, c('markers', 'locus')]
        pmap <- haplotypes.map[[name]][, c('markers', 'locus')]
        # anchors in pmap can't be assumed to be at the same loci as in dmap.
        pmap$locus <- .realign.loci(dmap, pmap)
        gm <- merge(dmap, pmap, all=T, sort=F)
        .debug.out(nrow(dmap), "markers in dmap,", nrow(pmap), "markers in pmap,", nrow(gm), "markers in merge.\n")
        .debug.out(length(unique(gm$markers)), "actual unique markers.\n")
        # there seems to still be a problem with abusive sorting. Damn R.
        gm <- gm[sort.int(gm$locus, na.last=NA, index.return=T)$ix, ]
        # now recompute dist.from and dist.to
        gm$dist.from = c(Inf, diff(gm$locus))
        gm$dist.to = c(diff(gm$locus), Inf)
        # publish merges map
        attr(gm, 'chrom.name') <- name
        attr(gm, 'doubles') <- .find.doubles(gm)
        #merged.maps[[name]] <- list(map=gm, doubles=.find.doubles(gm))
        merged.maps[[name]] <- gm
    }
    merged.maps
}


