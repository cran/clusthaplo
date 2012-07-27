#source('draw.chromosome.R')

.use.deprecated.draw.chromosomes <- F

#plot.clusthaplo.haplotypes <- function(clusters, ord=1:ncol(clusters), ...) {
plot.clusthaplo.haplotypes <- function(x, ...) {
    sup <- list(...)
    if(is.null(sup$ord)) {
        sup$ord <- 1:ncol(x)
    }
    if(.use.deprecated.draw.chromosomes) {
        if(is.null(sup$names)) {
            sup$names <- attr(x, 'ind.names')[sup$ord]
        }
        if(is.null(sup$loci)) {
            sup$loci <- attr(x, 'loci')
        }
        sup$cl <- .build.cliq.struc(.cliques(x))
        do.call(.draw.chromosome.OBSOLETE, sup)
        title(xlab="Map (cM)")
    } else {
        sup$tc <- x
        sup$step.size <- max(diff(attr(x, 'loci')))
        do.call(.draw.chrom, sup)
    }
}


.generic.hist.hack <- function(loci, counts, ylab, tick.interval=25, ...) {
    locus.max <- loci[length(loci)]
    count.max <- max(counts)

    int.dens <- rep(1 / length(loci), length(loci))

    pos.dif <- diff(loci)

    hi <- list(breaks=loci,
               counts=counts,
               intensities=int.dens,
               #density=int.dens,
               density=counts,  # so R doesn't whine with wrong surfaces. We just don't use freq=T
               mids=pos.dif / 2 + loci[-length(loci)],
               xname=as.character(loci),
               equidist=FALSE)
    class(hi) <- 'histogram'

    plot.new()
    sup.args <- list(...)
    pw.args <- list(xlim=c(0, locus.max), ylim=c(0, count.max), mar=c(0, 0, 0, 0), oma=c(0, 0, 0, 0), xlab="Map (cM)", ylab=ylab)
    for(n in names(sup.args)) {
        pw.args[[n]] <- sup.args[[n]]
    }
    do.call(plot.window, pw.args)
    .ticks <- unique(c(seq(0, locus.max, by=tick.interval),
                       locus.max))
    .labels <- as.character(round(.ticks, 1))
    axis(1,  # bottom
         labels=.labels,
         line=.1,
         at=.ticks)
    .ticks <- 0:count.max
    axis(2, at=.ticks)
    #lines(hi, col="gray40", border="gray40", freq=T)
    lines(hi, col="gray40", border="gray40")
    title(xlab="Map (cM)")
    #lines(loci, counts, t='h', lwd=c(1, diff(loci)))
}


plot.clusthaplo.density <- function(x, ...) {
    #dens <- dens$density
    #dens[is.na(dens)] <- 0
    .generic.hist.hack(x$locus, x$density, "Window density", ...)
}


plot.clusthaplo.alleles <- function(x, ...) {
    .generic.hist.hack(x$locus, x$count, "Number of ancestral alleles", ...)
}


.plot.clusters <- function(loci, clustering, ymin, ymax, ...) {
    polyx <- vector()
    polyx[(1:length(loci)) * 2 - 1] <- loci
    polyx[(2:length(loci)) * 2 - 2] <- loci[-1]
    polyx[2 * length(loci)] <- loci[length(loci)]
    polyx <- c(polyx, rev(polyx))
    polyy <- rep(ymin, length(polyx))
    polyy[(1:length(clustering)) * 2 - 1] <- ymin + clustering * (ymax - ymin)
    polyy[(1:length(clustering)) * 2] <- ymin + clustering * (ymax - ymin)
    polygon(polyx, polyy, ...)
}


plot.clusthaplo.pair.similarity <- function(x, ...) {
    ymax <- max(1, max(x$similarity), na.rm=T)
    plot.new()
    sup.args <- list(...)
    p.args <- list(x=x$locus, y=x$similarity, type="l",
                   xlab="Map (cM)", ylab="Similarity",
                   ylim=c(0, ymax), mar=c(0, 0, 0, 0), oma=c(0, 0, 0, 0))
    for(n in names(sup.args)) {
        p.args[[n]] <- sup.args[[n]]
    }
    do.call(plot, p.args)
    #plot(x$locus, x$similarity, type='l', xlab="Position", ylab="Similarity", ylim=c(0, ymax), ...)
    #polyx <- vector()
    #polyx[(1:length(x$locus)) * 2 - 1] <- x$locus
    #polyx[(2:length(x$locus)) * 2 - 2] <- x$locus[-1]
    #polyx[2 * length(x$locus)] <- x$locus[length(x$locus)]
    #polyx <- c(polyx, rev(polyx))
    #polyy <- rep(0, length(polyx))
    #polyy[(1:length(x$clustering)) * 2 - 1] <- x$clustering
    #polyy[(1:length(x$clustering)) * 2] <- x$clustering
    #polygon(polyx, polyy, col="#00FF0040", border="#00000000")
    .plot.clusters(x$locus, x$clustering, ymin=0, ymax=1, col="#00FF0040", border="#00000000")
    if(!is.null(attr(x, 'threshold'))) {
        abline(h=attr(x, 'threshold'), col='red')
    }
    #title()
}

